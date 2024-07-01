# Install

``` bash
git clone
cd hashtable_for_bqf
g++ main.cpp additionnal_methods.cpp -o main -I<unordered_dense_include_dir>
```

# Build

``` bash
./main build <input_counted_kmers_filepath> <output_binary_index_filepath>
```

# Query

``` bash
./main query <input_binary_index_filepath> <input_queries_filepath> <output_results_filepath> 
```
