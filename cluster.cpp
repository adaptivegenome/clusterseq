/*********************************************************************
 *
 * cluster.cpp
 * Part of the Clusterseq tool
 *
 * Author: Lee C. Baker, VBI
 * Last modified: 11 July 2012
 *
 * See https://github.com/adaptivegenome/clusterseq for more information.
 *
 *********************************************************************
 *
 * This file is released under the Virginia Tech Non-Commercial 
 * Purpose License. A copy of this license has been provided as
 * LICENSE.txt.
 *
 *********************************************************************//*
 *
 */
#include <cassert>
#include <cstdio>
#include <cstring>
#include <climits>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <errno.h>
#include <algorithm>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#else
#warning "Compiling without OpenMP. Enabling OpenMP will increase performance (add -fopenmp to your command line)."
#endif

using namespace std;

int cluster_distance(const string & s1, const string & s2, const int max_diff = INT_MAX) {
    size_t length = min(s1.size(), s2.size());
    int score = 0;
    
    for( size_t i = 0; i < length && score <= max_diff; i++)
        if(s1[i] != s2[i] && s1[i] != 'N' && s2[i] != 'N')
            score++;
    
    return score;
}

class step2_compare{
    map<string,int> & m;
public:
    step2_compare( map<string,int> & m)
    : m(m)
    {}

    bool operator() (const string & i,const string & j) const { return m[i] < m[j];}
};

void merge_clusters(vector<string> filenames, int min_count_for_filters, int cluster_edit_distance_threshold);

//find the most common letter for each position in the sequence. All must be the same length.
string getConsensus(vector<string> sequences)
{
    if(sequences.size() == 1)
        return sequences.front();
    size_t length = sequences.begin()->size();

    for(vector<string>::const_iterator i = sequences.begin(); i != sequences.end(); i++) {
        if(i->size() != length)
            cerr << "Sequences are clustered together with different lengths. Aborting." << endl;
    }
    
    static int count[256];  //store counters for each letter in here
    string output(length,'?');

    for(size_t ix = 0; ix < length; ix++) {
        memset(count, 0, sizeof(count));

        // accumulate counts for each letter
        // The consensus shouldn't include an N unless there is no other option, so we
        // increment scores by 2, and then limit N to a score of one.
        for(vector<string>::const_iterator i = sequences.begin(); i != sequences.end(); i++)
            count[int((*i)[ix])] += 2;

        count['N'] = min(1,count['N']);

        //find the index in the array of the most common letter.
        int * max_el = max_element(count, &count[sizeof(count) / sizeof(*count)]);
        char result_c = distance(count, max_el);

        output[ix] = result_c;
    }

    return output;
}

map<string, set<string> > readTagFile(string filename = "lala")
{
    map<string, set<string> > tags;
    
    ifstream file(filename.c_str());
    
    string last_name;
    while (!file.eof() && !file.fail()) {
        string line, tag;
        getline(file, line);
        stringstream line_str(line);
        
        if(file.eof() || file.fail())
            break;

        if(line[0] != '\t' && line[0] != ' ')
            line_str >> last_name;
        line_str >> tag;

        if(line_str.fail())
            break;
        
        tags[last_name].insert(tag);
    }
    
    cerr << tags.size() << " records loaded from tag file " << filename << "." << endl;
    
    return tags;
}

double frand() {
    return (double) rand() / (double) RAND_MAX;
}

int main(int cargs, char ** vargs)
{
    srand(time(NULL));
    if(cargs < 2) {
        cerr << "Filters and clusters data in a supplied FASTQ file." << endl;
        cerr << "Data is expected to be in the following format:" << endl;
        cerr << "  [tag][start_marker][data][end_marker]" << endl;
        cerr << endl;
        cerr << "Usage:" << endl;
        cerr << "    cluster [-a min_quality max_n_allowed num_diff_allowed] [-f sample_fraction] [-k known_barcodes] [-cf min_count_for_filters] [-cd cluster_edit_distance_threshold] FASTQ_name_no_ext [tag_file_name] [start_marker] [end_marker]" << endl;
        cerr << "\nFirst stage:\n" << endl;
        cerr << "    min_quality    : Minimum allowed quality. Bases with lower quality become 'N' (Default 30)" << endl;
        cerr << "    max_n_allowed  : Highest number of allowed 'N's per sequence- others discarded (Default 3)" << endl;
        cerr << "    num_diff_alwd  : Number of differences allowed between sequences clustered together (Default 3)" << endl;
        cerr << "    sample_fraction: Fraction of the input to keep, from 0 to 1.0 as a decimal (default 1.0)" << endl;
        cerr << "    known_barcode  : A file with known correct barcodes to cluster to (otherwise, cluster)" << endl;
        cerr << "    FASTQ_name     : Name of the input file (omit .fastq)" << endl;
        cerr << "    tag_file_name  : Name of the file listing tags for this input file" << endl;
        cerr << "    start_marker   : Sequence to look for at end of each read" << endl;
        cerr << "    end_marker     : Sequence to look for at beginning of read after tag" << endl;
        cerr << "\nSecond stage:\ngenerates merged_clusters.csv merged_clusters_filtered.csv merged_clusters_histogram.csv\n" << endl;
        cerr << "    cluster_edit_distance_threshold : Maximum allowable edit distance between tag clusters to be grouped (default 3)" << endl;
        cerr << "    min_count_for_filters           : Lines in merged_clusters_filtered.csv must have a sequence occurring at least this many times to be included in the file (default 1)" << endl;
        cerr << endl;
        cerr << "Options must be specified in the above order!" << endl;
        cerr << endl;
        return -1;
    }

    int carg_counter = 1;
    
    int min_quality = 63;
    int max_n_allowed = 3;
    int score_threshold = 3;

    if(0 == strcmp("-a", vargs[carg_counter])) {
        carg_counter++;
        min_quality = (int)strtol(vargs[carg_counter++], NULL, 10) + 33;

        if((0 == min_quality && errno == EINVAL) || min_quality < 33 || min_quality > 127) {
            cerr << "Error parsing minimum quality parameter." << endl;
            return -1;
        }
        
        max_n_allowed = (int)strtol(vargs[carg_counter++], NULL, 10);

        if((0 == max_n_allowed && errno == EINVAL) || max_n_allowed < 0) {
            cerr << "Error parsing max N parameter." << endl;
            return -1;
        }
        
        score_threshold = (int)strtol(vargs[carg_counter++], NULL, 10);

        if((0 == score_threshold && errno == EINVAL) || score_threshold < 0) {
            cerr << "Error parsing score threshold parameter." << endl;
            return -1;
        }
    }
    

    double keep_fraction = 1.;
    if(0 == strcmp("-f", vargs[carg_counter])) {
        carg_counter++;
        char * convert_out = NULL;
        keep_fraction = strtod(vargs[carg_counter++], &convert_out);
        carg_counter++;
    }
    
    vector<string> known_barcodes;
    if(0 == strcmp("-k", vargs[carg_counter])) {
        carg_counter++;
        
        ifstream file(vargs[carg_counter++]);

        while(true) {
            string line;
            getline(file, line);
            
            if(!file.fail())
                known_barcodes.push_back(line);
            else
                break;
        }
        
        cerr << "Read " << known_barcodes.size() << " known barcodes from file." << endl;
    }
    
    //stage two parameters:
    int min_count_for_filters = 1;
    int cluster_edit_distance_threshold = 3;
    if(0 == strcmp("-cf", vargs[carg_counter])) {
        carg_counter++;
        char * mc = NULL;
        min_count_for_filters = strtol(vargs[carg_counter++], &mc,10);
    }
    
    if(0 == strcmp("-cd", vargs[carg_counter])) {
        carg_counter++;
        char * mc = NULL;
        cluster_edit_distance_threshold = strtol(vargs[carg_counter++], &mc,10);
    }

    cerr << "Keeping sequences with quality of at least " << (min_quality -33) << " (ASCII " << (int) min_quality <<"='" << (char)min_quality << "') and max " << max_n_allowed << " 'N's." << endl;

    const string fastq_name(vargs[carg_counter++]);
    istream * infile = &cin;
    ifstream infile_actual;
    infile_actual.open(string(fastq_name + ".fastq").c_str());
    infile = &infile_actual;

    string tag_file_name("lala");
    if(carg_counter < cargs)
        tag_file_name = vargs[carg_counter++];

    map<string, set<string> > tag_file = readTagFile(tag_file_name);

    if(infile_actual.fail()) {
        cerr << "Error opening input file " << fastq_name << ".fastq" << endl;
        return -1;
    }

    if(tag_file.count(fastq_name) == 0) {
        cerr << "No entries for tag " << fastq_name << " in tag file " << tag_file_name << ". Aborting." << endl;
        return -1;
    }
    
    const set<string> & tags = tag_file[fastq_name];
    size_t tag_length = tags.begin()->size();
    for(std::set<std::string>::const_iterator tag_it = tags.begin(); tag_it != tags.end(); tag_it++) {
        if(tags.begin()->size() != tag_it->size()) {
            std::cerr << "Tags '" << *(tags.begin()) << "' and '" << *tag_it << "' have different lengths. Aborting." << std::endl;
            return -1;
        }
    }

    cerr << tag_file.size() << " FASTQ records loaded from tag file " << tag_file_name << "; using " << tags.size() << " tags for this FASTQ." << endl;
    string begin_marker("GGCGCGCC");
    
    if(carg_counter < cargs)
        begin_marker = vargs[carg_counter++];

    string end_marker;
    
    if(carg_counter < cargs)
        end_marker = vargs[carg_counter++];
    else {    
        switch(tag_length) {
            case 2:
                end_marker = "GCGGCC";
                break;
            case 4:
                end_marker = "GCGG";
                break;
            default:
                cerr << "Invalid tag length " << tag_length << ". Aborting." << endl;
                break;
        }
    }

    //read in data
    size_t seen_sequences = 0;
    size_t discarded_sequences = 0;
    size_t discarded_sequences_due_to_random = 0;
    size_t discarded_sequences_due_to_tags = 0;
    size_t discarded_sequences_due_to_markers = 0;
    size_t kept_sequences = 0;
    string seq, qual;
    const size_t data_start = begin_marker.size() + tag_length;
    
    //initialize tag output files
    map<string, ofstream *> tag_output_files;
    map<string, vector<string> *> sequences_map;
    for(set<string>::iterator i = tags.begin(); i != tags.end(); i++) {
        const string filename = fastq_name + "." + *i + ".txt";
        tag_output_files[*i] = new ofstream(filename.c_str());
        sequences_map.insert(pair<string, vector<string> *>(*i, new vector<string>()));
    }

    //process input file
    while(true) {
        
        seq.clear();
        qual.clear();
        //script one
        char name_line_start_sequence, name_line_start_quality;

        infile->get(name_line_start_sequence);
        infile->ignore(INT_MAX, '\n');
        getline(*infile, seq);
        infile->get(name_line_start_quality);
        infile->ignore(INT_MAX, '\n');
        getline(*infile, qual);

        if(infile->fail() || infile->eof())
            break;

        if(name_line_start_sequence != '@' || name_line_start_quality != '+') {
            cerr << "Skipping 4-line sequence in FASTQ file, bad format" << endl;
            continue;
        }
        
        seen_sequences++;
        const string tag = seq.substr(0,tag_length);

	if(frand() > keep_fraction) {
	    discarded_sequences_due_to_random++;
	    continue;
	}

        if(0 == tags.count(tag)) {
            discarded_sequences_due_to_tags++;
            continue;
        }

        //check for begin/end strings at correct positions
        if(0 != strncmp(begin_marker.c_str(), &(seq.c_str()[tag_length]), begin_marker.size())
            || 0 != strncmp(end_marker.c_str(), &(seq.c_str()[seq.size() - end_marker.size()]), end_marker.size())) {
                discarded_sequences_due_to_markers++;
            continue;
        }
        
        const size_t data_size = seq.size() - begin_marker.size() - end_marker.size() - tag.size();
        seq = seq.substr(data_start, data_size);
        qual = qual.substr(data_start, data_size);

        if(seq.size() < 2 || qual.size() < 2)
            break;
        
        *tag_output_files[tag] << seq << "\t" << qual << endl;

        //script 2
        if(seq.size() != qual.size()) {
            cerr << "WARNING: Skipping sequence because sequence length != quality length" << endl;
            continue;
        }

        //replace low quality reads with Ns
        for(size_t i = 0; i < seq.size(); i++) {
            if(qual[i] < min_quality)
                seq[i] = 'N';
        }

        if(count(seq.begin(), seq.end(), 'N') <= max_n_allowed) {
            kept_sequences++;
            sequences_map[tag]->push_back(seq);
        } else
            discarded_sequences++;
    }

    for(set<string>::iterator i = tags.begin(); i != tags.end(); i++)
        delete tag_output_files[*i];

    cerr << setw(8) << seen_sequences << " sequences read." << endl; 
    cerr << setw(8) << discarded_sequences_due_to_random << " (" << 100. * discarded_sequences_due_to_random / seen_sequences << "%) discarded due to random threshold(" << keep_fraction << ")." << endl;    
    cerr << setw(8) << discarded_sequences_due_to_tags << " (" << 100. * discarded_sequences_due_to_tags / seen_sequences << "%) discarded for invalid tags." << endl;    
    cerr << setw(8) << discarded_sequences_due_to_markers << " (" << 100. * discarded_sequences_due_to_markers / seen_sequences << "%) discarded for invalid begin or end markers." << endl; 
    cerr << setw(8) << discarded_sequences << " (" << 100. * discarded_sequences / seen_sequences << "%) discarded for too many Ns" << endl;
    cerr << setw(8) << kept_sequences << " (" << 100. * kept_sequences / seen_sequences << "%) kept." << endl;

    if(kept_sequences == 0) {
        cerr << "Didn't find any sequences that look like:" << endl;
        for(set<string>::iterator tag_it = tags.begin(); tag_it != tags.end(); tag_it++)
            cerr << "  " << *tag_it << begin_marker << "<data>" << end_marker << endl;
        
        cerr << "Check your begin and end markers, tags, and data." << endl;
        
    }
    
    vector <string> merge_cluster_filenames;

    for(set<string>::iterator tag_it = tags.begin(); tag_it != tags.end(); tag_it++) {
        vector<string> & sequences = *sequences_map[*tag_it];
        string tag = *tag_it;

        if(sequences.empty())
            continue;

        cerr << "Processing tag " << tag << "(" << sequences.size() << " sequences):" << endl;

        //step one - count unique sequences
        map<string, int> unique_sequence_counts;
        for(vector<string>::iterator i = sequences.begin(); i != sequences.end(); i++) {
            map<string, int>::iterator mi = unique_sequence_counts.find(*i);
            
            if(mi == unique_sequence_counts.end())
                unique_sequence_counts[*i] = 1;
            else
                mi->second++;    
        }

        cerr << setw(8) << unique_sequence_counts.size() << " unique sequences." << endl;
        
        //step two - sort sequences array by count
        vector<string> sorted_sequences;
        sorted_sequences.reserve(unique_sequence_counts.size());
        for(map<string, int>::const_iterator i = unique_sequence_counts.begin(); i != unique_sequence_counts.end(); i++)
            sorted_sequences.push_back(i->first);

        sort(sorted_sequences.begin(), sorted_sequences.end(), step2_compare(unique_sequence_counts));
        
        //step three - perform clustering
        multimap<string, string> cluster_members;
        vector<string> & cluster_centers = known_barcodes.empty() ? sorted_sequences : known_barcodes;

#ifdef _OPENMP
#pragma omp parallel shared(sorted_sequences,score_threshold,cluster_members) default(none)
        {
            multimap<string, string> cluster_members_local;
#pragma omp for
            for(int i = 0; i < (int)sorted_sequences.size(); i++) {
                const string & sequence = sorted_sequences[i];
                for(vector<string>::const_reverse_iterator j = cluster_centers.rbegin(); j != cluster_centers.rend(); j++) {
                    if(cluster_distance(sequence, *j, score_threshold) <= score_threshold) {
                        cluster_members_local.insert(pair<string, string>(*j,sequence));
                        break;
                    }
                }
            }
#pragma omp critical
            {
                cluster_members.insert(cluster_members_local.begin(), cluster_members_local.end());
            }
        }
#pragma omp barrier

#else
        for(size_t i = 0; i < sorted_sequences.size(); i++) {
            const string & sequence = sorted_sequences[i];
            for(vector<string>::const_reverse_iterator j = sorted_sequences.rbegin(); j != sorted_sequences.rend(); j++) {
                if(cluster_distance(sequence, *j, score_threshold) <= score_threshold) {
                    cluster_members.insert(pair<string, string>(*j,sequence));
                    break;
                }
            }
        }
#endif

        //step four - count cluster members
        //step five - build list of clusters (combined)
        map<string, int> cluster_count;
        
        for (multimap<string,string>::iterator j = cluster_members.begin(); j!=cluster_members.end(); j++) {
            const string & cluster_name = j->first;
            map<string, int>::iterator ci = cluster_count.find(cluster_name);
            if(ci == cluster_count.end())
                cluster_count[cluster_name] = unique_sequence_counts[j->second];
            else
                ci->second += unique_sequence_counts[j->second];
        }
        cerr << setw(8) << cluster_count.size() << " clusters." << endl;

        //step six- build consensus list. Build list of sequences to go in, then call getConsensus
        map<string, string> consensuses;    //maps cluster to consensus
        for(map<string,int>::const_iterator i = cluster_count.begin(); i != cluster_count.end(); i++) {
            vector<string> consensus_in;
            const string & cluster_name = i->first;

            for (multimap<string,string>::iterator j = cluster_members.equal_range(cluster_name).first; j!=cluster_members.equal_range(cluster_name).second; j++)
                for(int ctr = 0; ctr < unique_sequence_counts[j->second]; ctr++)
                    consensus_in.push_back(j->second);

            consensuses[cluster_name] = getConsensus(consensus_in);
        }
        
        merge_cluster_filenames.push_back(fastq_name + "." + tag + "_clusters.csv");
        ofstream outfile(merge_cluster_filenames.back().c_str());
        
        //output consensuses
        for(map<string, string>::const_iterator i = consensuses.begin(); i != consensuses.end(); i++)
            outfile << i->second << "," << cluster_count[i->first] << endl;
        cerr << setw(8) << consensuses.size() << " sequences written." << endl;
    }
    
    merge_clusters(merge_cluster_filenames, min_count_for_filters, cluster_edit_distance_threshold);

    return 0;
}

class value_sorter : less<pair<string, int> > {
public:
    bool operator() (const pair<string, int> & a, const pair <string,int> & b) {
        return a < b;
    }
};

string cluster_name(const string & a, const string & b) {
    assert(a.size() == b.size());
    string n = a;
    
    for(size_t i = 0; i < a.size(); i++)
        if(a[i] != b[i])
            n[i] = 'N';
    return n;
}

void merge_clusters(vector<string> filenames, int min_count_for_filters, int cluster_edit_distance_threshold) {
    vector<map<string, int> > file_data(filenames.size());
    set<string> keys;
    int ct = 0;
    double histogram_bin_growth_factor = 1.4;
    
    //load data from files
    for(vector<string>::const_iterator filename = filenames.begin(); filename != filenames.end(); filename++) {
        ifstream file(filename->c_str());
        
        string line;
        while(true) {
            string name;
            int count;
            getline(file, line);
            
            istringstream ss( line );
            getline( ss, name, ',' );
            ss >> count;
            
            if(!file.good() || ss.fail())
                break;
            
            file_data[ct][name] = count;
            file_data[ct].insert(pair<string, int>(name, count));
            keys.insert(name);
        }
        ct++;
    }
    
    //identify keys to be remapped
    map<string,int> counts;
    vector<pair<string,int> > sorted_keys;
    for(set<string>::const_iterator key = keys.begin(); key != keys.end(); key++) {
        counts[*key] = 0;
        for(vector<map<string,int> >::iterator data = file_data.begin(); data != file_data.end(); data++) {
            if(data->count(*key))
                counts[*key] += (*data)[*key];
        }
        sorted_keys.push_back(std::pair<string, int>(*key, counts[*key]));
    }
    
    map<string, string> remapped_keys;
    int attached = 0;
    
    value_sorter vs;
    sort(sorted_keys.begin(), sorted_keys.end(), vs);
    
    for(vector<pair<string,int> >::const_reverse_iterator i1 = sorted_keys.rbegin(); i1 != sorted_keys.rend(); i1++) {
        for(vector<pair<string,int> > ::const_iterator i2 = sorted_keys.begin(); i2 != sorted_keys.end(); i2++) {
            if(i1->first == i2->first)
                break;
            if(cluster_distance(i1->first, i2->first, cluster_edit_distance_threshold+2) <= cluster_edit_distance_threshold) {
                remapped_keys[i2->first] = i1->first;
                if(!remapped_keys.count(i1->first))
                    remapped_keys[i1->first] = i1->first;
                //cerr << "Attaching " << i2->first << " to cluster " << i1->first << endl;
                attached += 1;
            }
        }
    }
    
    //remove keys mapped to a cluster node that doesn't exist any more
    set<string> keys_seen_once, keys_seen_twice, keys_difference;
    
    for(map<string, string>::const_iterator i = remapped_keys.begin(); i != remapped_keys.end(); i++) {
        if(!keys_seen_once.count(i->second)) {
            keys_seen_once.insert(i->second);
            continue;
        }
        if(!keys_seen_twice.count(i->second))
            keys_seen_twice.insert(i->second);
    }
    
    
    set_difference(keys_seen_once.begin(), keys_seen_once.end(), keys_seen_twice.begin(), keys_seen_twice.end(), inserter(keys_difference, keys_difference.end()));
    int removed = 0;
    vector<string> to_remove;
    for(map<string, string>::const_iterator i = remapped_keys.begin(); i != remapped_keys.end(); i++) {
        if(!keys_difference.count(i->second)) {
            //cerr << "Removing orphan cluster mapping " << i->first << " to " << i->second << endl;
            //remapped_keys.erase(i->first);
            to_remove.push_back(i->first);
            removed++;
        }
    }
    for(vector<string>::const_iterator i = to_remove.begin(); i != to_remove.end(); i++)
        remapped_keys.erase(*i);
    
    cerr << "Clustering stage 2: Attached " << attached << " sequences to clusters" << endl;
    
    // generate new names for cluster centers
    map<string, string> remapped_names;
    for(map<string, string>::const_iterator i = remapped_keys.begin(); i != remapped_keys.end(); i++) {
        if(remapped_names.count(i->second))
            remapped_names[i->second] = cluster_name(i->first, remapped_names[i->second]);
        else
            remapped_names[i->second] = i->second;
    }
    
    for(map<string, string>::const_iterator i = remapped_keys.begin(); i != remapped_keys.end(); i++)
        cerr << "Renaming cluster " << i->first << " as " << i->second << endl;
    
    //now perform merging of files
    vector<map<string, int> > new_file_data;
    
    for(vector<map<string,int> >::iterator data = file_data.begin(); data != file_data.end(); data++) {
        new_file_data.push_back(map<string, int>());
        map<string,int> & new_data = new_file_data.back();
        
        for(set<string>::const_iterator key = keys.begin(); key != keys.end(); key++) {
            if(data->count(*key)) {
                if(remapped_keys.count(*key)) {
                    string & cluster_key = remapped_keys[*key];
                    string & cluster_name = remapped_names[cluster_key];
                    if(new_data.count(cluster_name))
                        new_data[cluster_name] += (*data)[*key];
                    else
                        new_data[cluster_name] = (*data)[*key];
                } else
                    new_data[*key] = (*data)[*key];
            }
        }
    }
    
    file_data = new_file_data;
    
    //generate the new keys list after clustering
    keys.clear();
    for(vector<map<string,int> >::iterator data = file_data.begin(); data != file_data.end(); data++) {
        for(map<string, int>::const_iterator i = data->begin(); i != data->end(); i++)
            keys.insert(i->first);
    }
    
    //generate merged counts output file
    {
        stringstream header;
        header << "sequence";
        for(vector<string>::const_iterator i = filenames.begin(); i != filenames.end(); i++)
            header << "," << *i;
        ofstream outfile("merged_clusters.csv");
        ofstream outfile_filtered("merged_clusters_filtered.csv");
        outfile << header.str();
        outfile_filtered << header.str();
        
        for(set<string>::const_iterator key = keys.begin(); key != keys.end(); key++) {
            bool keep = false;
            stringstream line;
            line << *key;
            for(vector<map<string,int> >::iterator data = file_data.begin(); data != file_data.end(); data++) {
                if(data->count(*key)) {
                    line << "," << (*data)[*key];
                    
                    if((*data)[*key] > min_count_for_filters)
                        keep = true;
                } else {
                    line << ",0";
                }
            }
            outfile << line.str() << "\n";

            if(keep)
                outfile_filtered << line.str() << "\n";
        }
        outfile_filtered.close();
    }
    
    //generate histogram
    {
        //generate list of sizes (counts for each key)
        vector<int> sizes;
        for(set<string>::const_iterator key = keys.begin(); key != keys.end(); key++) {
            int ct = 0;
            
            for(vector<map<string,int> >::iterator data = file_data.begin(); data != file_data.end(); data++) {
                if(data->count(*key))
                    ct += (*data)[*key];
            }
            sizes.push_back(ct);
        }
        
        sort(sizes.begin(), sizes.end());
        
        ofstream outfile;
        outfile.open("merged_clusters_histogram.csv");
        double bin_size(1.0);
        
        //generate bins for histogram
        vector<int> bins;
        bins.push_back(1);
        while(bin_size < sizes[sizes.size() - 1]) {
            double new_bin_size = bin_size * histogram_bin_growth_factor;
            if(int(bin_size) != int(new_bin_size))
                bins.push_back(int(bin_size));
            bin_size = new_bin_size;
        }
        bins.push_back(int(bin_size));
        
        int bin_ct = 0, s = 0;
        
        //generate actual histogram
        for(size_t i = 0; i < sizes.size(); i++) {
            while(sizes[i] > bins[bin_ct + 1]) {
                if(bins[bin_ct] != bins[bin_ct + 1]) {
                    outfile << bins[bin_ct] << "-" << bins[bin_ct + 1] << "," << s << "\n";
                    s = 0;
                }
                
                bin_ct += 1;
            }
            s += sizes[i];
        }
        outfile << bins[bin_ct] << "-" << bins[bin_ct + 1] << "," << s << endl;
    }
}
