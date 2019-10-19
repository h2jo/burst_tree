# burst_tree_decomposition

Author: Hang-Hyun Jo (h2jo23@gmail.com)

Upload date: 19 October 2019

Description: 
    Codes developed and used for the following paper: Hang-Hyun Jo, Takayuki Hiraoka, and Mikko Kivela, Burst-tree decomposition of time series reveals the structure of temporal correlations [arXiv:1907.13556] 


Folder description

1. /header: header files

2. /burst_tree: basic analysis using burst-tree decomposition method

    main.c: main c file
    
    main.x: main execution file
    
    run.sh: shell script for compiling and running the program
    
    timings_example.txt: the first 1000 events for the most active Wikipedia editor
    
    /sample_result: the results files using "timings_example.txt"

3. /MRRM: microcanonical randomized reference models

    main.c: main c file
    
    main.x: main execution file
    
    run.sh: shell script for compiling and running the program
    
    timings_example.txt: the first 1000 events for the most active Wikipedia editor
    
    tree_example.txt: the tree structure derived from "timings_example.txt"

4. /model: kernel-based modeling
    
    main.c: main c file
    
    main.x: main execution file
    
    run.sh: shell script for compiling and running the program
