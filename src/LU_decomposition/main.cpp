#include <omp.h>
#include <string>
#include <filesystem>
#include <vector>
#include <cstdio> // Для printf

#include "tests.cpp"
#include "lu_decomposition.cpp"

// Путь к каталогу для сохранения результатов
std::filesystem::path PATH_TO_RESULTS = "results/M1_8_core";

int main() {
    // Получение максимального количества потоков
    int max_threads = omp_get_max_threads();
    printf("[INFO] Max threads: %d\n", max_threads);
    
    int num_threads = 1;
    std::vector<int> n_vec = {512, 768, 1024, 1280, 1536, 1792, 2048};
    
    int block_size = 128;

    // Запуск всех тестов для LU-разложения с использованием блочного метода
    make_all_tests_LU(
        num_threads, 
        PATH_TO_RESULTS / ("n" + std::to_string(num_threads) + "b" + std::to_string(block_size) + "lu_block.txt"), 
        n_vec, 
        make_parallel_block_LU_decomposition, 
        block_size
    );

    return 0;
}
