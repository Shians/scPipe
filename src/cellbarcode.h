#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <memory>
#include <asio/thread_pool.hpp>
#include <asio/post.hpp>
#include <Rcpp.h>
#include "utils.h"

#ifndef CELLBARCODE_H
#define CELLBARCODE_H
// a class that stores cellular barcode annotation and
// find close barcode for a given sequence
class Barcode
{
public:
    std::unordered_map<std::string, std::string> barcode_dict;
    std::vector<std::string> cellid_list;
    std::vector<std::string> barcode_list;

    // if annotation is given
    void read_anno(std::string fn);
    void set_threads(int nthreads);

    std::unordered_map<std::string, std::string> get_count_file_path(std::string out_dir);

    std::string get_closest_match(const std::string &bc_seq, int max_mismatch);

    friend std::ostream& operator<< (std::ostream& out, const Barcode& obj);
private:
    asio::thread_pool *cell_bc_t_pool_;
    int n_threads_;
};

#endif
