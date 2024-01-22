# Solving the Max-Cut problem with Neutral Atoms

The repo is organized as follows:

- Main notebooks explaining all the functions and answering the questions and points from the coding challenge:
  - QAOA_maxcut_neutral_atoms
  - VQAA_maxcut_neutral_atoms
- Bonus work testing the VQAA and QAOA with different graph instances:
  - QAOA_tests
  - VQAA_tests
  - encoding_tests
- The bonus works with 3 different modules that are inside the folder src:
  - mapping (functions to obtain the coordinates for the encoding into the atomic register)
  - qubo_tools (functions to create the QUBO matrix for the Max-Cut problem)
  - quantum_solver_qaoa
  - quantum_solver_vqaa
- I used Pre-commit to mantain certain standard in the code. Pre-commit applies a uniform style on everyting and will check for basic problems in the code. 
The rest of the file are for the functioning of the Pre-commit process.