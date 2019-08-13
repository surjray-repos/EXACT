%% example of how to use Ancestree

% input file
input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input';
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% parameters
n_mutations = size(F_from_SampleData,1);
T_samples = size(F_from_SampleData,2);
disp(n_mutations);
disp(T_samples);

scale = scaling;
alpha = 0.3; % if alpha is big, lots of things will be clustered together
beta = 0.8; % to choose a larger beta we need more samples , i.e. larger T_samples
gamma = 0.01; % small gamma means larger confidence on the data
wrapper_working_directory = [pwd, '/distribution/'];
ancestree_executable_path = [pwd, '/distribution/build/ancestree'];

% run AncesTree
ancestree_output = ancestree_wrapper(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);

% get output
%since AncesTree might have multiple solutions, we shall check how many
%solutions are inferred and reference them by sol_id
sol_id = 1; %Considering the first opimal tree
Tree_Matrix_T = ancestree_output{3}{sol_id};
Mutant_Frequencies_M = ancestree_output{2}{sol_id};

