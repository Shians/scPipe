// demultiplexing
#include "cellbarcode.h"

using std::string;
using std::unordered_map;
using std::unique_ptr;
using namespace Rcpp;

// read annotation from files, the file should have two columns
// first column is cell id and second column is barcode.
// protocols with two barcodes are merged as one.
void Barcode::read_anno(string fn)
{
    check_file_exists(fn);
    std::ifstream infile(fn);
    string line;
    std::getline(infile, line); // skip header

    char sep;
    if (line.find(",") != string::npos)
    {
        sep = ',';
    }
    else if (line.find("\t") != string::npos)
    {
        sep = '\t';
    }
    else
    {
        Rcpp::stop("the annotation file should be comma or tab separated");
    }

    while(std::getline(infile, line))
    {
        std::stringstream linestream(line);
        if (line.size() < 3){
            continue;
        }
        string cell_id;
        string barcode;

        std::getline(linestream, cell_id, sep);
        std::getline(linestream, barcode, sep);

        barcode_dict[barcode] = cell_id;

        if (std::find(cellid_list.begin(), cellid_list.end(), cell_id) == cellid_list.end())
        {
            cellid_list.push_back(cell_id);
        }
        
        barcode_list.push_back(barcode);
    }
}

void Barcode::set_threads(int nthreads) {
    n_threads_ = std::max(nthreads, 1);
    cell_bc_t_pool_ = new asio::thread_pool(n_threads_);
}

unordered_map<string, string> Barcode::get_count_file_path(string out_dir)
{
    string csv_fmt = ".csv";
    unordered_map<string, string> out_fn_dict;
    for(const auto& n : cellid_list) 
    {
        out_fn_dict[n] = join_path(out_dir, n+csv_fmt);
    }
    return out_fn_dict;
}

string Barcode::get_closest_match(const string &bc_seq, int max_mismatch)
{
    if (barcode_dict.find(bc_seq) != barcode_dict.end())
    {
        return bc_seq;
    }
    int sml1st = std::numeric_limits<int>::max();
    int sml2ed = std::numeric_limits<int>::max();
    std::vector<int> hamming_dists(barcode_list.size());
    string closest_match;

    int tasks_per_thread = 0;
    int chunk_start = 0;
    if (n_threads_ > 1) {
        tasks_per_thread = barcode_list.size() / (n_threads_ - 1);
        const int CHUNK_SIZE = std::max(256, tasks_per_thread);
        int chunk_end = CHUNK_SIZE;

        while (chunk_end < barcode_list.size()) {
            asio::dispatch(*cell_bc_t_pool_, 
                [&, chunk_start, chunk_end] () {
                    for (int i = chunk_start; i < chunk_end; i++) {
                        hamming_dists[i] = hamming_distance(barcode_list[i], bc_seq);
                    }
                }
            );
            chunk_start += CHUNK_SIZE;
            chunk_end += CHUNK_SIZE;
        }
    }

    for (int i = chunk_start; i < barcode_list.size(); i++)
    {
        hamming_dists[i] = hamming_distance(barcode_list[i], bc_seq);
    }

    cell_bc_t_pool_->join();

    for (int i = 0; i < hamming_dists.size(); i++)
    {
        string const &bc = barcode_list[i];
        int dist = hamming_dists[i];
        if (dist <= max_mismatch)
        {
            if (dist < sml1st)
            {
                sml1st = dist;
                closest_match = bc;
            }
            else if (sml1st <= dist && dist < sml2ed)
            {
                sml2ed = dist;
            }
        }
    }

    if (sml1st < sml2ed)
    {
        return closest_match;
    }
    else
    {
        return string();
    }

}

std::ostream& operator<< (std::ostream& out, const Barcode& obj)
{
    for( const auto& n : obj.barcode_dict ) 
    {
        out << "Barcode:[" << n.first << "] Cell Id:[" << n.second << "]\n";
    }
    return out;
}
