# Public Transport Route Planning: When Connectivity and Commuting Demand Both Matter
This repo holds the code and dataset for the reproducibility of CT-Bus.

## Running
1. Function test: we can run a specific configuration on NYC, with k=15, sn=5000, Tn=1, and w=0.5.

./fairbus_run.m

2. Submodularity test: we can reproduce Figure 2.

./test_submodular.m

3. Sensibility test: We can reproduce Figures 6 and 8.

./experiments.m


## Dataset

Chicago dataset is in folder "chi_data".

NYC dataset is in folder "nyc_data".

Each folder has four files, storing the potential edges, existing edges, node location, and a random matrix for connectivity computation.

*** We also share the dataset that are not needed for running but used to generate the above dataset.

Trajectory dataset can be downloaded here: https://drive.google.com/open?id=1NFUH-YHOXDoSAYZqJPwVkJK8E-H56TIw

## Results
After running successfully, our code will generate two files: "chi_busroutes_new.txt" is ids of stops, and another one "chi_route_stops.txt" generates a set of locations.

## Visualization

Users can visualize based on files "chi_fairbus_k501.txt" and "nyc_fairbus_k501.txt" using MapV (https://mapv.baidu.com/) library by opening file "show_fairbus_nyc" and "show_fairbus_chi" using WebStorm (https://www.jetbrains.com/webstorm/).
