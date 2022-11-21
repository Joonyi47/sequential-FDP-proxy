# Sequential-FDP-proxy
Sequential robust FDP optimization with CNN and LSTM proxy models.

## References
- [Sequential field development plan through robust optimization coupling with CNN and LSTM-based proxy models](https://doi.org/10.1016/j.petrol.2021.109887)
- [Robust optimization of the locations and types of multiple wells using CNN based proxy models](https://doi.org/10.1016/j.petrol.2020.107424)
- [Efficient and robust optimization for well patterns using a PSO algorithm with a CNN-based proxy model](https://doi.org/10.1016/j.petrol.2021.109088) 

##
You can do the optimizations with the following steps. (Folder names represent each step)

## 1.	Well placement

### 1.1.	Proxy modeling (CNN-TOF)

1) Move to the directory for the target field
2) Run ‘main_for_data_sampling.m’
    1. The permeability models are stored in the ‘PERM.mat’
    2. Indexes for the representative models are stored in the ‘PERM_selected_idx.mat’
    3. The number of sample data depends on the value in the [Np] variable
    4. The operation conditions depend on the values in the [wset] variable
    5. Set file name for saving sampled data at the end of the code
   * you can download example [sample data](https://drive.google.com/file/d/1K1OI9K-0PSkzUMoxnXKcVvDtKsHR7Pyi/view?usp=share_link)
   
3) Run ‘main_for_model_training.m’ (recommended to run each section sequentially)
    1. In Initial conditions section, set [varname] variable same as the file name
    2. In Initial conditions section, set [Nall] variable same as the number of sample data
    3. In Input data section, the input data is incorporated as different ways according to the [case_sensitivity] variable
    4. In CNN model section, set structure and hyperparameters of the CNN model and train the model
    5. In later sections, visualize regression plots for the training, validation, and test data and save the CNN model
4) Copy ‘trainedNet.mat’ in the defined directory to the "\1.2. Optimization\target field\data"

### 1.2.	Optimization

1) Move to the directory for the target field
2) Run ‘main_for_PSOwCNNup.m’
    1. In the first section, load the CNN model in ‘trainedNet.mat’
    2. For several runs, label each run with the [casename] variable
    3. The generations for the retraining process depend on the value in the [rtgen] variable
3) Copy ‘pos_opt”casename”.mat’ to the "\2.1. Proxy modeling (LSTM)\target field\data"



## 2.	Well operation

### 2.1.	Proxy modeling (LSTM)

1) Move to the directory for the target field
2) Run ‘main_for_data_sampling.m’
    1. In the first section, load the well placement opt. result in ‘pos_opt[casename].mat’
    2. The number of sample data depends on the value in the [Np] variable
    3. The ranges of operation conditions depend on the values in the [wset] variable 
3) Set file name for saving sampled data at the end of the code
4) Run ‘main_for_model_training.m’ (recommended to run each section sequentially)
    1.	In Data preprocessing section, data structure is transformed to be used in the LSTM model
    2.	In Initial conditions section, set [varname] variable same as the file name
    3.	In Initial conditions section, set [Nall] variable same as the number of sample data
    4.	In LSTM model section, set structure and hyperparameters of the LSTM model and train the model
    5.	In later sections, visualize regression plots for the training, validation, and test data and save the LSTM model
5) Copy ‘trainedNet.mat’ in the defined directory to the "\2.2. Optimization\target field\data"

### 2.2.	Optimization

1) Move to the directory for the target field
2) Run ‘main_for_PSOwLSTMup.m’
    1. In the first section, load the LSTM model in ‘trainedNet.mat’
    2. In the first section, lead the well placement opt. result in ‘pos_opt[casenam].mat’
    3. For several runs, label each run with the [casename] variable
    4. The generations for the retraining process depend on the value in the [rtgen] variable
3) Get the final solution in ‘pos_opt[casename].mat’
