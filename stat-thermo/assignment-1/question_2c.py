import numpy as np

# Define the sample size (we can change this to whatever)
sample_size = 10000

# Generate random values between -1 and 3
# To do this, use the formula: a + (b - a) * np.random.random(size)
# where a is the lower bound and b is the upper bound
x_uniform = -1 + (3 - (-1)) * np.random.random(sample_size)

# Open the file to write the generated values
with open('x_uniform.dat', 'w') as f:
    for x in x_uniform:
        f.write(f"{x}\n")  # Writing each random number to the file

print("Random numbers between -1 and 3 have been generated and written to 'x_uniform.dat'.")