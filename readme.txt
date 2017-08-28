Usage: Call the program with the following command line parameters: <mode> <filename>. where mode is one of the following: 0 : normal operation 1 : control 2 : decimal to binary convertor"
Control is used to check if a given solution (a list of correction terms) satisfies the non-completeness, correctness and uniformity constraints

main.cpp is the main file
quadratic.cpp is an unfinished version that includes quadratic correction terms
/Input contains the input in binary form, this is the form used by the program in mode 0
/DecimalInput contains the input in decimal form, this can be converted to binary form by using mode 2
/output contains the solutions found by the program

example usage:
Convert to binary form:
$  ./uniformity.out 2 input_5_28.in
run:
$  ./uniformity.out 0 input_5_28.in
