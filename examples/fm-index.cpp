#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace sdsl;
using namespace std;

bool sortbyfirst(const tuple<int,string>& a,  
               const tuple<int,string>& b) 
{ 
    return (get<0>(a) < get<0>(b)); 
} 


void replace_char(vector<char>* extracted_vec, vector<char> sigma, char c){
    for (unsigned i = 0; i < sigma.size(); i++) {
        if (!(sigma.at(i) == c)){
            (*extracted_vec).push_back(sigma.at(i));
        }
	}
    
}


int main(int argc, char** argv)
{
    if (argc <  2) {
        cout << "Usage " << argv[0] << " text_file [max_locations] [post_context] [pre_context]" << endl;
        cout << "    This program constructs a very compact FM-index" << endl;
        cout << "    which supports count, locate, and extract queries." << endl;
        cout << "    text_file      Original text file." << endl;
        cout << "    max_locations  Maximal number of location to report." <<endl;
        cout << "    post_context   Maximal length of the reported post-context." << endl;
        cout << "    pre_context    Maximal length of the pre-context." << endl;
        return 1;
    }
    string error_num = argv[2];
    size_t max_locations = 5;
    size_t post_context = 10;
    size_t pre_context = 10;
    if (argc >= 4) {
        max_locations = atoi(argv[3]);
    }
    if (argc >= 5) {
        post_context = atoi(argv[4]);
    }
    if (argc >= 6) {
        pre_context = atoi(argv[5]);
    }
    string index_suffix = ".fm9";
    string index_file   = string(argv[1])+index_suffix;
    string rev_string_prefix = "rev_";
    string rev_file = rev_string_prefix + string(argv[1]);
    string rev_index_file   = string(rev_file)+index_suffix;
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_index;
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> rev_fm_index;

    std::ifstream ifs(argv[1]);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                       (std::istreambuf_iterator<char>()));
    
    int input_length = content.length();
    reverse(content.begin(), content.end());
    content = content.substr(1, content.length()-1);
    std::ofstream outfile (rev_file);
    outfile << content << std::endl;
    outfile.close();

    if (!load_from_file(fm_index, index_file)) {
        ifstream in(argv[1]);
        if (!in) {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index "<<index_file<< " located. Building index now." << endl;
        construct(fm_index, argv[1], 1); // generate index
        store_to_file(fm_index, index_file); // save it
       
    }
    if (!load_from_file(rev_fm_index, rev_index_file)) {
        ifstream in(rev_file);
        if (!in) {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index "<<rev_index_file<< " located. Building index now." << endl;
        construct(rev_fm_index, rev_file, 1); // generate index
        store_to_file(rev_fm_index, rev_index_file); // save it
       
    }

    cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << endl;
    cout << "Input search terms and press Ctrl-D to exit." << endl;
    string prompt = "\e[0;32m>\e[0m ";
    cout << prompt;
    string query;
    while (getline(cin, query)) {
        vector<char> sigma;
        int sig_size = fm_index.sigma;
        for (int i = 1 ; i < sig_size ; i++){
            sigma.push_back(fm_index.comp2char[i]);
        }
        int q_len = query.length();
        char q_copy_temp[q_len];
        strcpy(q_copy_temp,query.c_str());
        string q_copy(q_copy_temp);
        if (error_num.compare("0") == 0){
            size_t m  = query.size();
            size_t occs = sdsl::count(fm_index, query.begin(), query.end());
            cout << "# of occurrences: " << occs << endl;
            if (occs > 0) {
                cout << "Location and context of first occurrences: " << endl;
                auto locations = locate(fm_index, query.begin(), query.begin()+m);
                sort(locations.begin(), locations.end());
                for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i) {
                    cout << setw(8) << locations[i] << ": ";
                    if (pre_extract > locations[i]) {
                        pre_extract = locations[i];
                    }
                    if (locations[i]+m+ post_extract > fm_index.size()) {
                        post_extract = fm_index.size()-locations[i]-m;
                    }
                    auto s   = extract(fm_index, locations[i]-pre_extract, locations[i]+m+ post_extract-1);
                    string pre = s.substr(0, pre_extract);
                    s = s.substr(pre_extract);
                    if (pre.find_last_of('\n') != string::npos) {
                        pre = pre.substr(pre.find_last_of('\n')+1);
                    }
                    cout << pre;
                    cout << "\e[1;31m";
                    cout << s.substr(0, m);
                    cout << "\e[0m";
                    string context = s.substr(m);
                    cout << context.substr(0, context.find_first_of('\n')) << endl;
                }
            }
        }
        else if (error_num.compare("1") == 0) {
            int num_of_occs = 0;
            vector<int> found;
            vector<tuple<int, string>> first_occ;
            size_t m  = query.size();
            int middle = m/2;
            
            std::string s_half;
            if (m % 2 == 0){
                s_half = query.substr (m/2,m/2);
            }
            else{
                s_half = query.substr (m/2,m/2 + 1);
            }
            


            for (int i = 0 ; i < middle ; i++){
                string acc1 = query.substr(middle - i, q_len - middle + i);
                auto sh_address = lex_interval(fm_index,query.begin() + (middle -i), query.end());
                uint64_t lb_i = sh_address[0];
                uint64_t rb_i = sh_address[1];
                char e = query[(middle - i -1)];
                vector<char> * sh_extracted = new vector<char>(0);
                replace_char(sh_extracted, sigma, e);

                int ext_size = sh_extracted->size();
                
                for (int j = 0 ; j < ext_size ; j++){
                    vector<int> curr_occ;
                    uint64_t lb_j = lb_i;
                    uint64_t rb_j = rb_i;
                    char temp_char = (*sh_extracted).at(j);
                    char* curr_char = &temp_char;
                    string acc2 = temp_char + acc1;
                    backward_search(fm_index, lb_j, rb_j, (char)*curr_char, lb_j, rb_j);
                    bool cont = true;
                    if (lb_j > rb_j){
                        cont = false;
                    }
                    
                    for (int k = middle - i - 2 ; k >= 0 && cont  ; k--){
                        
                        char temp_c = query[k];
                        char* c = &temp_c;
                        acc2 = temp_c + acc2;
                        backward_search(fm_index, lb_j, rb_j, (char)*c, lb_j, rb_j);
                        if (lb_j > rb_j){
                        cont = false;
                        }
                    }
                    if (cont){
                        if (max_locations == 0){
                            num_of_occs = num_of_occs + (rb_j -lb_j +1);
                        }
                        else{
                            for (unsigned l = lb_j ; l <= rb_j ; l++){
                                curr_occ.push_back(fm_index[l]);
                                found.push_back(fm_index[l]);
                            }
                            for (unsigned k =0 ; k < curr_occ.size() ; k++){
                                first_occ.push_back(make_tuple(curr_occ.at(k), acc2));
                            }
                            curr_occ.clear();
                        }
                    }
                    
                    

                }
            delete(sh_extracted);
            }
            
            // f_half = first half of the input string query
            reverse(query.begin(),query.end());
            std::string f_half;
            
            int middle2;
            if (m % 2 == 0){
                middle2 = m/2;
            }
            else{
                middle2 = m/2 + 1;
            }
            if (m % 2 == 0){
                f_half = query.substr (middle,m/2);
            }
            else{
                f_half = query.substr (middle,m/2);
            }
            
            for (int i = 0 ; i < middle2 ; i++){
                string acc1 = query.substr(middle - i + 1, q_len - middle + i);
                auto sh_address = lex_interval(rev_fm_index,query.begin() + (middle2 -i), query.end());
                uint64_t lb_i = sh_address[0];
                uint64_t rb_i = sh_address[1];
                char e = query[(middle2 - i -1)];
                vector<char> * sh_extracted = new vector<char>(0);
                replace_char(sh_extracted, sigma, e);

                int ext_size = sh_extracted->size();
                
                for (int j = 0 ; j < ext_size ; j++){
                    vector<int> curr_occ;
                    uint64_t lb_j = lb_i;
                    uint64_t rb_j = rb_i;
                    char temp_char = (*sh_extracted)[j];
                    char* curr_char = &temp_char;
                    string acc2 = temp_char + acc1;
                    backward_search(rev_fm_index, lb_j, rb_j, (char)*curr_char, lb_j, rb_j);
                    bool cont = true;
                    if (lb_j > rb_j){
                        cont = false;
                    }
                    
                    for (int k = middle2 - i - 2 ; k >= 0 && cont  ; k--){
                        char temp_c = query[k];
                        char* c = &temp_c;
                        acc2 = temp_c + acc2;
                        backward_search(rev_fm_index, lb_j, rb_j, (char)*c, lb_j, rb_j);
                        if (lb_j > rb_j){
                        cont = false;
                        }
                    }
                    
                    if (cont){
                        if (max_locations == 0){
                            num_of_occs = num_of_occs + (rb_j -lb_j +1);
                        }
                        else {
                            for (unsigned l = lb_j ; l <= rb_j ; l++){
                                int res = rev_fm_index[l];
                                found.push_back(input_length - res -query.length() -1);
                                curr_occ.push_back(input_length - res -query.length() -1);
                            }
                            reverse(acc2.begin(),acc2.end());
                            for (unsigned k =0 ; k < curr_occ.size() ; k++){
                                first_occ.push_back(make_tuple(curr_occ.at(k), acc2));
                            }
                            curr_occ.clear();
                        }
                    }
                }
            delete(sh_extracted);
            }
            int found_size = found.size();
            int size = max(found_size,num_of_occs);
            cout << q_copy << " : " << size << endl;
        
            int all_occ_num = min(max_locations, first_occ.size());
            if (all_occ_num > 0){
                cout << "Location and context of first occurrences: "<< endl;
                sort(first_occ.begin(), first_occ.end(), sortbyfirst);
            }
            for (int i = 0; i < all_occ_num ; i = i + 1){
                cout << "  " << get<0>(first_occ.at(i)) << ": " << get<1>(first_occ.at(i)) << endl;
            }

        }
        else{
            int num_of_occs = 0;
            vector<tuple<int,string>> all_occ;
            int s1 = floor(q_len/3);
            int s2 = q_len -s1;
            // CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A
            vector<int> found_a;
            vector<tuple<int,string>> first_occ_a;
            for (int t = s2 -1  ; t > 0 ; t--){ // t is e1 index
                string acc1 = query.substr(t + 1, q_len - t - 1);
                auto sh_address = lex_interval(fm_index,query.begin() + t + 1 , query.end());
                uint64_t lb_i = sh_address[0];
                uint64_t rb_i = sh_address[1];
                char e1 = query[t];
                vector<char> * sh_extracted = new vector<char>(0);
                replace_char(sh_extracted, sigma, e1);
                int ext_size1 = sh_extracted->size();       
                for (int j = 0 ; j < ext_size1 ; j++){
                    uint64_t lb_j = lb_i;
                    uint64_t rb_j = rb_i;
                    char temp_char = (*sh_extracted)[j];
                    char* curr_char = &temp_char;
                    backward_search(fm_index, lb_j, rb_j, (char)*curr_char, lb_j, rb_j);
                    bool cont = true;
                    if (lb_j > rb_j){
                        cont = false;
                    }
                    for (int k = t-1 ; k >= 0 && cont  ; k--){
                        string acc2 = temp_char + acc1;
                        auto sh_address1 = lex_interval(fm_index,acc2.begin() , acc2.end());
                        uint64_t lb_k = sh_address1[0];
                        uint64_t rb_k = sh_address1[1];
                        for (int l = t-1 ; l > k ; l-- ){
                            char temp_c = query[l];
                            char* c = &temp_c;
                            acc2 = temp_c + acc2;
                            backward_search(fm_index, lb_k, rb_k, (char)*c, lb_k, rb_k);
                            if (lb_k > rb_k){
                            cont = false;
                            }
                        }
                        char e2 = query[k]; //// k is e2 index
                        vector<char> * sh_extracted = new vector<char>(0);
                        replace_char(sh_extracted, sigma, e2);
                        int ext_size2 = sh_extracted->size();

                        for (int p = 0 ; p < ext_size2 ; p++){
                            cont = true;
                            uint64_t lb_p = lb_k;
                            uint64_t rb_p = rb_k;
                            char temp_char = (*sh_extracted)[p];
                            char* curr_char = &temp_char;
                            string acc3 = temp_char + acc2;
                            backward_search(fm_index, lb_p, rb_p, (char)*curr_char, lb_p, rb_p);
                            if (lb_p > rb_p){
                                cont = false;
                            }
                            for (int r = k - 1 ; r >= 0 && cont  ; r--){
                                
                                char temp_c = query[r];
                                char* c = &temp_c;
                                acc3 = temp_c + acc3;
                                
                                backward_search(fm_index, lb_p, rb_p, (char)*c, lb_p, rb_p);
                                if (lb_p > rb_p){
                                cont = false;
                                }
                            }
                            vector<int> curr_occ;
                            if (cont){
                                if (max_locations == 0){
                                   num_of_occs = num_of_occs + (rb_p -lb_p +1);
                                }
                                else{
                                    for (unsigned z = lb_p ; z <= rb_p ; z++){
                                        curr_occ.push_back(fm_index[z]);
                                        found_a.push_back(fm_index[z]);
                                    }
                                    for (unsigned k =0 ; k < curr_occ.size() ; k++){
                                        first_occ_a.push_back(make_tuple(curr_occ.at(k), acc3));
                                        all_occ.push_back(make_tuple(curr_occ.at(k), acc3));
                                    }
                                    curr_occ.clear();
                                }
                            }
                            cont = true;   
                        }
                        delete(sh_extracted);
                    }
                    
                }
                delete(sh_extracted);
            }
            

            // CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B
            vector<int> found_b;
            vector<tuple<int,string>> first_occ_b;
            reverse(query.begin(), query.end());
            int temp_s2 = s1; 
            for (int t = temp_s2 -1  ; t > 0 ; t--){ // t is e1 index
                string acc1 = query.substr(t + 1, q_len - t - 1);
                auto sh_address = lex_interval(rev_fm_index,query.begin() + t + 1 , query.end());
                uint64_t lb_i = sh_address[0];
                uint64_t rb_i = sh_address[1];
                char e1 = query[t];
                vector<char> * sh_extracted = new vector<char>(0);
                replace_char(sh_extracted, sigma, e1);
                int ext_size1 = sh_extracted->size();       
                for (int j = 0 ; j < ext_size1 ; j++){
                    uint64_t lb_j = lb_i;
                    uint64_t rb_j = rb_i;
                    char temp_char = (*sh_extracted)[j];
                    char* curr_char = &temp_char;
                    backward_search(rev_fm_index, lb_j, rb_j, (char)*curr_char, lb_j, rb_j);
                    bool cont = true;
                    if (lb_j > rb_j){
                        cont = false;
                    }
                    for (int k = t-1 ; k >= 0 && cont  ; k--){
                        string acc2 = temp_char + acc1;
                        auto sh_address1 = lex_interval(rev_fm_index,acc2.begin(), acc2.end());
                        uint64_t lb_k = sh_address1[0];
                        uint64_t rb_k = sh_address1[1];
                        for (int l = t-1 ; l > k ; l-- ){
                            char temp_c = query[l];
                            char* c = &temp_c;
                            acc2 = temp_c + acc2;
                            backward_search(rev_fm_index, lb_k, rb_k, (char)*c, lb_k, rb_k);
                            if (lb_k > rb_k){
                            cont = false;
                            }
                        }
                        char e2 = query[k]; //// k is e2 index
                        vector<char> * sh_extracted = new vector<char>(0);
                        replace_char(sh_extracted, sigma, e2);
                        int ext_size2 = sh_extracted->size();
                        for (int p = 0 ; p < ext_size2 ; p++){
                            cont = true;
                            uint64_t lb_p = lb_k;
                            uint64_t rb_p = rb_k;
                            char temp_char = (*sh_extracted)[p];
                            char* curr_char = &temp_char;
                            string acc3 = temp_char + acc2;

                            backward_search(rev_fm_index, lb_p, rb_p, (char)*curr_char, lb_p, rb_p);
                            if (lb_p > rb_p){
                                cont = false;
                            }
                            for (int r = k - 1 ; r >= 0 && cont  ; r--){
                                
                                char temp_c = query[r];
                                char* c = &temp_c;
                                acc3 = temp_c + acc3;
                                
                                backward_search(rev_fm_index, lb_p, rb_p, (char)*c, lb_p, rb_p);
                                if (lb_p > rb_p){
                                    cont = false;
                                }
                            }
                            
                            reverse(acc3.begin(),acc3.end());
                            vector<int> curr_occ;
                            if (cont){
                                if (max_locations == 0){
                                    num_of_occs = num_of_occs + (rb_p -lb_p +1);
                                }
                                else
                                {
                                    for (unsigned z = lb_p ; z <= rb_p ; z++){
                                        int res = rev_fm_index[z];
                                        found_b.push_back(input_length - res -query.length() -1);
                                        curr_occ.push_back(input_length - res -query.length() -1);
                                    }
                                    for (unsigned k =0 ; k < curr_occ.size() ; k++){
                                        first_occ_b.push_back(make_tuple(curr_occ.at(k), acc3));
                                        all_occ.push_back(make_tuple(curr_occ.at(k), acc3));
                                    }
                                    curr_occ.clear();
                                }
                            
                                
                            }
                            cont = true;   
                        }
                        delete(sh_extracted);
                    }
                    
                }
                delete(sh_extracted);
            }


            // CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C
            string passer;
            uint64_t lb_passer;
            uint64_t rb_passer;
            int third;
            int two_thirds;
            if (q_len % 3 == 1){
                third = q_len/3;
                two_thirds = third*2 + 2;
            }
            else if (q_len % 3 == 2)
            {
                third = q_len/3 + 1;
                two_thirds = third*2 + 1;
            }
            else
            {
                third = q_len/3;
                two_thirds = third * 2 + 1;
            }
            third = q_len - third;
            two_thirds = q_len - two_thirds + 1;
            vector<int> found_c;
            vector<tuple<int,string>> first_occ_c;
            for (int t = third - 1  ; t > two_thirds - 1 ; t--){ // t is e1 index
                string acc1 = query.substr(t + 1, q_len - t - 1);
                auto sh_address = lex_interval(rev_fm_index,query.begin() + t + 1, query.end());
                uint64_t lb_i = sh_address[0];
                uint64_t rb_i = sh_address[1];
                char e1 = query[t];
                vector<char> * sh_extracted = new vector<char>(0);
                replace_char(sh_extracted, sigma, e1);
                int ext_size1 = sh_extracted->size();       
                for (int j = 0 ; j < ext_size1 ; j++){
                    uint64_t lb_j = lb_i;
                    uint64_t rb_j = rb_i;
                    char temp_char = (*sh_extracted)[j];
                    char* curr_char = &temp_char;
                    string acc2 = temp_char + acc1;
                    backward_search(rev_fm_index, lb_j, rb_j, (char)*curr_char, lb_j, rb_j);
                    bool cont = true;
                    if (lb_j > rb_j){
                        cont = false;
                    }
                    for (int r = t - 1 ; r > two_thirds - 1 && cont  ; r--){
                                
                        char temp_c = query[r];
                        char* c = &temp_c;
                        acc2 = temp_c + acc2;
                        
                        backward_search(rev_fm_index, lb_j, rb_j, (char)*c, lb_j, rb_j);

                        if (lb_j > rb_j){
                        cont = false;
                        }
                    }
                    passer = acc2;
                    lb_passer = lb_j;
                    rb_passer = rb_j;
                    
                    cont = true;
                    for (int z = two_thirds - 1 ; z >= 0 ; z--){ // t is e1 index
                        string acc3 = passer;
                        uint64_t lb_z = lb_passer;
                        uint64_t rb_z = rb_passer;
                        for (int w = two_thirds -1 ; w > z ; w--){
                            char temp_c = query[w];
                            char* c = &temp_c;
                            acc3 = temp_c + acc3;
                                
                            backward_search(rev_fm_index, lb_z, rb_z, (char)*c, lb_z, rb_z);
                            if (lb_z > rb_z){
                                cont = false;
                            }
                        }
                        
                        char e1 = query[z];
                        vector<char> * sh_extracted = new vector<char>(0);
                        replace_char(sh_extracted, sigma, e1);
                        int ext_size1 = sh_extracted->size();
                        for (int n = 0 ; n < ext_size1 ; n++){
                            uint64_t lb_n = lb_z;
                            uint64_t rb_n = rb_z;
                            char temp_char = (*sh_extracted)[n];
                            char* curr_char = &temp_char;
                            string acc4 = temp_char + acc3;
                            backward_search(rev_fm_index, lb_n, rb_n, (char)*curr_char, lb_n, rb_n);
                            bool cont = true;
                            if (lb_n > rb_n){
                                cont = false;
                            }
                            for (int r = z - 1 ; r >= 0 && cont  ; r--){
                                        
                                char temp_c = query[r];
                                char* c = &temp_c;
                                acc4 = temp_c + acc4;
                                
                                backward_search(rev_fm_index, lb_n, rb_n, (char)*c, lb_n, rb_n);
                                if (lb_n > rb_n){
                                cont = false;
                                }
                            }
                            vector<int> curr_occ;
                            if (cont){
                                if (max_locations == 0){
                                    num_of_occs = num_of_occs + (rb_n -lb_n +1);
                                }
                                else {
                                    reverse(acc4.begin(), acc4.end());
                                    for (unsigned z = lb_n ; z <= rb_n ; z++){
                                        int res = rev_fm_index[z];
                                        found_c.push_back(input_length - res -query.length() -1);
                                        curr_occ.push_back(input_length - res -query.length() -1);
                                    }
                                    for (unsigned k =0 ; k < curr_occ.size() ; k++){
                                        first_occ_c.push_back(make_tuple(curr_occ.at(k), acc4));
                                        all_occ.push_back(make_tuple(curr_occ.at(k), acc4));
                                    }
                                    curr_occ.clear();
                                }
                            }
                            cont = true; 
                        }
                    delete(sh_extracted);
                    }
                }
            delete(sh_extracted);  
            }
            
            

            // CASE D CASE D CASE D CASE D CASE D CASE D CASE D CASE D CASE D CASE D CASE D CASE D CASE D 
            reverse(query.begin(), query.end());
            third = q_len - third;
            two_thirds = q_len - two_thirds;
            int dist = two_thirds - third - 1;
            int two_thirds_rev  = q_len - two_thirds;
            uint64_t sp_lb;
            uint64_t sp_rb;
            vector<int> found_d;
            vector<tuple<int,string>> first_occ_d;
            
            
            for (int t = third - 1 ; t >= 0 ; t--){ // t is e1 index
                dist += 1;
                string acc1 = query.substr(t + 1, dist);
                auto sh_address1 = lex_interval(fm_index, acc1.begin(), acc1.end());
                uint64_t initial_lb = sh_address1[0];
                uint64_t initial_rb = sh_address1[1];
                
              

                char e1 = query[t];
                vector<char> * sh_extracted = new vector<char>(0);
                replace_char(sh_extracted, sigma, e1);
                int ext_size1 = sh_extracted->size();
                
                for (int j = 0 ; j < ext_size1 ; j++){
                    bool cont = true;
                    uint64_t lb_j = initial_lb;
                    uint64_t rb_j = initial_rb;
                    char temp_char = (*sh_extracted)[j];
                    char* curr_char = &temp_char;
                    string acc2 = temp_char + acc1;
                    backward_search(fm_index, lb_j, rb_j, (char)*curr_char, lb_j, rb_j);
                    if (lb_j > rb_j){
                        cont = false;
                    }
                    
                    for (int r = t - 1 ; r >= 0 && cont  ; r--){
                                
                        char temp_c = query[r];
                        char* c = &temp_c;
                        acc2 = temp_c + acc2;
                        
                        backward_search(fm_index, lb_j, rb_j, (char)*c, lb_j, rb_j);

                        if (lb_j > rb_j){
                        cont = false;
                        }
                    }
                    reverse(acc2.begin(), acc2.end());
                    auto sh_address2 = lex_interval(rev_fm_index, acc2.begin(), acc2.end());
                    uint64_t rlb = sh_address2[0];
                    uint64_t rrb = sh_address2[1];
                    reverse(acc2.begin(), acc2.end());
                    
                    passer = acc2;
                    rb_passer = rrb;
                    lb_passer = rlb;
                    sp_rb = rrb;
                    sp_lb = rlb;
                    
                    
                    for (int k = two_thirds_rev - 1 ; k >= 0 ; k--){ // t is e1 index
                        reverse(query.begin(),query.end());
                        if (k != two_thirds_rev -1){
                            char add = query[k + 1];
                            char* c = &add;
                            passer = passer + add;
                            backward_search(rev_fm_index, sp_lb, sp_rb, (char)*c, sp_lb, sp_rb);
                            if (sp_lb > sp_rb){
                                cont = false;
                            }
                            rb_passer = sp_rb;
                            lb_passer = sp_lb;

                        }
                        string acc3 = passer;
                        char e1 = query[k];
                        vector<char> * sh_extracted = new vector<char>(0);
                        replace_char(sh_extracted, sigma, e1);
                        int ext_size1 = sh_extracted->size();
                        for (int w = 0 ; w < ext_size1 ; w++){
                            bool cont = true;
                            uint64_t lb_w = lb_passer;
                            uint64_t rb_w = rb_passer;
                            char temp_char = (*sh_extracted)[w];
                            char* curr_char = &temp_char;
                            string acc4 = acc3  + temp_char;
                            backward_search(rev_fm_index, lb_w, rb_w, (char)*curr_char, lb_w, rb_w);
                            if (lb_w > rb_w) {
                                cont = false;
                            }
                            for (int z = k - 1 ; z >= 0 && cont  ; z--){
                                
                                char temp_c = query[z];
                                char* c = &temp_c;
                                acc4 = acc4 + temp_c;
                                
                                backward_search(rev_fm_index, lb_w, rb_w, (char)*c, lb_w, rb_w);

                                if (lb_w > rb_w){
                                    cont = false;
                                }
                            }
                            vector<int> curr_occ;
                            if (cont){
                                if (max_locations == 0){
                                    num_of_occs = num_of_occs + (rb_w -lb_w +1);
                                }
                                else {
                                for (unsigned z = lb_w ; z <= rb_w ; z++){
                                    int res = rev_fm_index[z];
                                    found_d.push_back(input_length - res -query.length() -1);
                                    curr_occ.push_back(input_length - res -query.length() -1);
                                }
                                for (unsigned k =0 ; k < curr_occ.size() ; k++){
                                    first_occ_d.push_back(make_tuple(curr_occ.at(k), acc4));
                                    all_occ.push_back(make_tuple(curr_occ.at(k), acc4));
                                }
                                curr_occ.clear();
                                }
                            }

                         }
                        reverse(query.begin(),query.end());
                        delete(sh_extracted);
                    }

                } 
            delete(sh_extracted);  
            }

            int found_size = found_a.size() + found_b.size() + found_c.size() + found_d.size();
            int size = max(found_size,num_of_occs);
            cout << q_copy << " : " << size << endl;

            int all_occ_num = min(max_locations, all_occ.size());
            if (all_occ_num > 0){
                cout << "Location and context of first occurrences: "<< endl;
                sort(all_occ.begin(), all_occ.end(), sortbyfirst);
            }
            for (int i = 0; i < all_occ_num ; i = i + 1){
                cout << "  " << get<0>(all_occ.at(i)) << ": " << get<1>(all_occ.at(i)) << endl;
            }

        }
        cout << prompt;
    }
    cout << endl;

}

