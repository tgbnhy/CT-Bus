# EquaBus-code
This repo holds the code and dataset for reproducibility of EquaBus.

## Running
1. Function test

./fairbus_run.m

2. Sumdularity test

./test_submodular.m

3. Sensibility test

./experiments.m


## Dataset

Chicago dataset is in "chi_data".

NYC dataset is in "nyc_data".

Each folder has four files, storing the potential edges, existing edges, node location, and a random matrix for connectivity computation.

## Results
After running successfully, our code will generate two files: "chi_busroutes_new.m" is ids of stops, and another one "chi_route_stops" generate a set of locations.

## Visualization

Users can visualize based on files "chi_fairbus_k501.txt" and "nyc_fairbus_k501.txt" using MapV (https://mapv.baidu.com/) library by opening file "show_fairbus_nyc" and "show_fairbus_chi" using WebStorm (https://www.jetbrains.com/webstorm/).
