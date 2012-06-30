/*
 * ClusterSeq
 * (c) 2012 Virginia Bioinformatics Institute
 *
 * Author: Lee Baker, VBI
 * See https://github.com/adaptivegenome/clusterseq for more information.
 */
#include <cassert>
#include <cstdio>
#include <cstring>
#include <climits>
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

#ifdef USE_OPENMP
#include <omp.h>
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

        if(line[0] != '\t')
            line_str >> last_name;
        line_str >> tag;
        
        tags[last_name].insert(tag);
    }
    
    cerr << tags.size() << " records loaded from tag file " << filename << "." << endl;
    
    return tags;
}

int main(int cargs, char ** vargs)
{
    if(cargs != 5 && cargs != 6 && cargs != 7 && cargs != 8) {
        cerr << "Filters and clusters data in a supplied FASTQ file." << endl;
        cerr << "Data is expected to be in the following format:" << endl;
        cerr << "  [tag][start_marker][data][end_marker]" << endl;
        cerr << endl;
        cerr << "Usage:" << endl;
        cerr << "    cluster min_quality max_n_allowed num_diff_allowed FASTQ_name_no_ext [tag_file_name] [start_marker] [end_marker]" << endl;
        cerr << endl;
        cerr << "    min_quality  : Minimum allowed quality. Bases with lower quality become 'N'" << endl;
        cerr << "    max_n_allowed: Highest number of allowed 'N's per sequence- others discarded" << endl;
        cerr << "    num_diff_alwd: Number of differences allowed between sequences clustered together" << endl;
        cerr << "    FASTQ_name   : Name of the input file (omit .fastq)" << endl;
        cerr << "    tag_file_name: Name of the file listing tags for this input file" << endl;
        cerr << "    start_marker : Sequence to look for at end of each read" << endl;
        cerr << "    end_marker   : Sequence to look for at beginning of read after tag" << endl;
        return -1;
    }

    int min_quality;
    if(isdigit(vargs[1][0])) {
        min_quality = (int)strtol(vargs[1], NULL, 10);
        if(min_quality < 33 && min_quality > 0)
            min_quality += 33;
    } else
        min_quality = (int) vargs[1][0];

    if((0 == min_quality && errno == EINVAL) || min_quality < 33 || min_quality > 127) {
        cerr << "Error parsing minimum quality parameter." << endl;
        return -1;
    }
    
    int max_n_allowed = (int)strtol(vargs[2], NULL, 10);

    if((0 == max_n_allowed && errno == EINVAL) || max_n_allowed < 0) {
        cerr << "Error parsing max N parameter." << endl;
        return -1;
    }
    
    int score_threshold = (int)strtol(vargs[3], NULL, 10);

    if((0 == score_threshold && errno == EINVAL) || score_threshold < 0) {
        cerr << "Error parsing score threshold parameter." << endl;
        return -1;
    }
    
    string tag_file_name("lala");
    if(cargs >= 6)
        tag_file_name = vargs[5];

    cerr << "Keeping sequences with quality of at least " << (min_quality -33) << " (ASCII " << (int) min_quality <<"='" << (char)min_quality << "') and max " << max_n_allowed << " 'N's." << endl;

    map<string, set<string> > tag_file = readTagFile(tag_file_name);
    istream * infile = &cin;
    ifstream infile_actual;

    infile_actual.open(string(string(vargs[4]) + ".fastq").c_str());
    infile = &infile_actual;

    if(infile_actual.fail()) {
        cerr << "Error opening input file " << vargs[4] << ".fastq" << endl;
        return -1;
    }

    const string fastq_name(vargs[4]);
    if(tag_file.count(fastq_name) == 0){
        cerr << "No entries for tag " << fastq_name << " in tag file. Aborting." << endl;
        return -1;
    }
    
    const set<string> & tags = tag_file[vargs[4]];
    size_t tag_length = tags.begin()->size();
    string begin_marker("GGCGCGCC");
    
    if(cargs >= 7)
        begin_marker = vargs[6];

    string end_marker;
    
    if(cargs == 8)
        end_marker = vargs[7];
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

#ifdef USE_OPENMP
#pragma omp parallel shared(sorted_sequences,score_threshold,cluster_members) default(none)
        {
            multimap<string, string> cluster_members_local;
#pragma omp for
            for(int i = 0; i < (int)sorted_sequences.size(); i++) {
                const string & sequence = sorted_sequences[i];
                for(vector<string>::const_reverse_iterator j = sorted_sequences.rbegin(); j != sorted_sequences.rend(); j++) {
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
        
        ofstream outfile((fastq_name + "." + tag + "_clusters.csv").c_str());
        
        //output consensuses
        for(map<string, string>::const_iterator i = consensuses.begin(); i != consensuses.end(); i++)
            outfile << i->second << "," << cluster_count[i->first] << endl;
        cerr << setw(8) << consensuses.size() << " sequences written." << endl;
    }

    return 0;
}
