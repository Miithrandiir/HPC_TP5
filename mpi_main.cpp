#include <iostream>
#include <mpi.h>
#include <random>
#include <iomanip>

#define NB_POINTS 10000
#define REAL_PI 3.141592653589793238462643
#define PRECISION 6
#define MASTER_PROC 0

void generateNumber();

int processId;
int nbProcess;

uint64_t totalLoop = 0;
uint64_t totalInCircle = 0;

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm groupe = MPI_COMM_WORLD;
    MPI_Comm_size(groupe, &nbProcess);
    MPI_Comm_rank(groupe, &processId);

    generateNumber();


    MPI_Finalize();
}

void generateNumber() {


    std::random_device rd;
    std::mt19937_64 mt(rd() + processId);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    bool needContinue = true;
    int nbLoop = NB_POINTS / nbProcess;

    while (needContinue) {
        int nbDansCercle = 0;

        for (int i = 0; i < nbLoop; i++) {

            double x = dist(mt);
            double y = dist(mt);
            if (x * x + y * y < 1) {
                ++nbDansCercle;
            }

        }

        int data[2]{nbLoop, nbDansCercle};

        int recvData[2];
        MPI_Reduce(data, recvData, 2, MPI_INT, MPI_SUM, MASTER_PROC, MPI_COMM_WORLD);

        //SI proc Maitre alors on regarde les résultats ...
        if (processId == MASTER_PROC) {

            totalLoop += recvData[0];
            totalInCircle += recvData[1];

            std::cout << "[MASTER] => recv Loop: " << recvData[0] << " \t recv Circle : " << recvData[1] << std::endl;
            std::cout << "[MASTER] => Total Loop: " << totalLoop << " \t Total In Circle : " << totalInCircle << std::endl;

            long double PI = 4.0 * ((long double) totalInCircle / (long double) totalLoop);
            if (std::abs(PI - REAL_PI) < std::pow(10, -PRECISION)) {
                needContinue = !needContinue;
                std::cout << "Found a good PI = " << PI << std::endl;
            } else {
                std::cout << "Bad PI  = " << std::scientific << std::setprecision(PRECISION) << PI << std::endl;
            }
            //broadcast
        }

        MPI_Bcast(&needContinue, 1, MPI_CXX_BOOL, MASTER_PROC, MPI_COMM_WORLD);
    }
}