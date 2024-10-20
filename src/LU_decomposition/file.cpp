#include <iostream>
#include <fstream>
#include <string>

/**
 * Записывает время time в файл.
 */
void write_time_to_txt(std::ofstream& file, double time) {
    file << std::to_string(time) << std::endl;
}
