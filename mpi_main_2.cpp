#include <iostream>
#include <mpi.h>
#include <random>
#include <iomanip>

#define NB_POINTS 10000
#define REAL_PI 3.141592653589793238462643
#define PRECISION 7
#define STOP_NUMBER -42

void generateNumber();

int SERVER_PROC = 1;
int MASTER_PROC = 0;

int processId;
int nbProcess;

int totalLoop = 0;
int totalInCircle = 0;

MPI_Comm main_comm;

void serverFunc();

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm groupe = MPI_COMM_WORLD;
    MPI_Comm_size(groupe, &nbProcess);
    MPI_Comm_rank(groupe, &processId);

    bool isServer = false;
    bool isMain = false;

    int color = 1;
    if (processId == SERVER_PROC) {
        color = 0;
    }

    MPI_Comm_split(MPI_COMM_WORLD, color, processId, &main_comm);

    //std::cout << processId << std::endl;

    if(processId == SERVER_PROC) {
        serverFunc();
    } else {
        generateNumber();
    }


    MPI_Finalize();

    return EXIT_SUCCESS;
}

void serverFunc() {
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<long double> dist(-1.0, 1.0);
    while (true) {
        int N;
        MPI_Status status;

        //std::cout << "[SERVER] Waiting for a Msg" << std::endl;
        MPI_Recv(&N, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        //std::cout << "[SERVER] Receive a msg with data : " << N << std::endl;

        //Check if it's a stop number
        if (N == STOP_NUMBER && status.MPI_SOURCE == MASTER_PROC) {
            break;
        }


        long double res[N];
        for (int i = 0; i < N; ++i) {
            res[i] = dist(mt);
        }

        //std::cout << "[SERVER] Send data to client" << std::endl;
        MPI_Send(res, N, MPI_LONG_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
    }
}

void generateNumber() {

    bool needContinue = true;
    int nbLoop = NB_POINTS;
    std::cout << "[CHILD] Proc Id " << processId << std::endl;
    while (needContinue) {
        int nbDansCercle = 0;
        //Demande des informations + traitements


        //std::cout << "[CHILD] Asking for numbers" << std::endl;
        MPI_Send(&nbLoop, 1, MPI_INT, SERVER_PROC, 0, MPI_COMM_WORLD);
        long double recvRandomData[nbLoop];
        //std::cout << "[CHILD] Waiting for numbers" << std::endl;
        MPI_Recv(recvRandomData, nbLoop, MPI_LONG_DOUBLE, SERVER_PROC, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //std::cout << "[CHILD] Numbers get" << std::endl;

        for (int i = 0; i < nbLoop; i += 2) {
            long double x = recvRandomData[i];
            long double y = recvRandomData[i + 1];
            if (x * x + y * y < 1) {
                ++nbDansCercle;
            }

        }
        //NbLoop/2 car on utilise 2 nombre aléatoire par boucle
        int data[2]{nbLoop/2, nbDansCercle};

        int recvData[2];
        std::cout << "[CHILD] Reducing" << std::endl;
        int res = MPI_Reduce(data, recvData, 2, MPI_INT, MPI_SUM, MASTER_PROC, main_comm);

        if (res != MPI_SUCCESS) {
            std::cout << "[ERROR] Not success reduce !" << std::endl;
        }

        //SI proc Maitre alors on regarde les résultats ...
        if (processId == MASTER_PROC) {

            //std::cout << "[MASTER] => recv Loop: " << recvData[0] << " \t recv Circle : " << recvData[1] << std::endl;
            totalLoop += recvData[0];
            totalInCircle += recvData[1];

            //std::cout << "[MASTER] => Total Loop: " << totalLoop << " \t Total In Circle : " << totalInCircle << std::endl;

            long double PI = 4.0 * ((long double) totalInCircle / (long double) totalLoop);
            if (std::abs(PI - REAL_PI) < std::pow(10, -PRECISION)) {
                needContinue = !needContinue;
                std::cout << "Found a good PI = " << PI << std::endl;
                int i = STOP_NUMBER;
                //Notification du serveur que l'on peut arrêter
                MPI_Send(&i, 1, MPI_INT, SERVER_PROC, 0, MPI_COMM_WORLD);
            } else {
                std::cout << "Bad PI  = " << std::scientific << std::setprecision(PRECISION) << PI << std::endl;
            }
            //broadcast
        }

        MPI_Bcast(&needContinue, 1, MPI_CXX_BOOL, MASTER_PROC, main_comm);
    }
}