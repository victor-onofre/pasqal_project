# Solving the Max-Cut problem with Neutral Atoms

In this work, I solved the Max-cut problem using the Quantum Approximate Optimization Algorithm (QAOA) and the Variational Quantum Adiabatic Algorithm (VQAA) 
with neutral atoms. All simulations were done with [Pulser](https://pulser.readthedocs.io/en/stable/). The report of the project is [here](https://github.com/victor-onofre/pasqal_project/blob/main/Report_PASQAL_modeling_coding_challenge_Victor_Onofre.pdf).

The repo is organized as follows:

- Main notebooks explaining all the functions and answering the questions and points from the coding challenge:
  - QAOA_maxcut_neutral_atoms (If you want to you can run it in [Colab here](https://colab.research.google.com/drive/1JQ6xAjDC7PWPLdLgyw5JcEEAhTYJC_az?usp=sharing) )
  - VQAA_maxcut_neutral_atoms ([Colab link](https://colab.research.google.com/drive/1vtuOBDd5rM7Ii0by-HJhUyXuNzvopmdh?usp=sharing))
- Bonus notebooks testing the VQAA and QAOA with different graph instances:
  - QAOA_tests
  - VQAA_tests
  - encoding_tests
- The bonus notebooks works with 3 different modules that are inside the folder src:
  - mapping (functions to obtain the coordinates for the encoding into the atomic register)
  - qubo_tools (functions to create the QUBO matrix for the Max-Cut problem)
  - quantum_solver_qaoa
  - quantum_solver_vqaa
- I used Pre-commit to maintain certain standards in the code. Pre-commit applies a uniform style on everything and will check for basic problems in the code.  The rest of the files are for the functioning of the Pre-commit process.
