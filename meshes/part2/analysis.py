import sys
import matplotlib.pyplot as plt
import numpy as np

# Get the array from command-line arguments
time_points = sys.argv[1:]

# Now you can use the execution_times array in your Python script
for i in range(len(time_points)):
    time_points[i] = float(time_points[i])

execution_times_short = []
execution_times_long = []
execution_times_large = []

r1 = np.arange(1, 21, 2)

for i in range(0, len(time_points), 2):
    exe = round((time_points[i + 1] - time_points[i]), 4)
    if i >= ((len(time_points) / 3) * 2):
        execution_times_short.append(exe)
    elif i >= (len(time_points) / 3):
        execution_times_long.append(exe)
    else:
        execution_times_large.append(exe)

fig = plt.figure(figsize=(30, 15))
width1 = 1
print(execution_times_large)
print(execution_times_short)
print(execution_times_long)

# creating the bar plot
plt.bar(r1, execution_times_short, width=0.5, color="blue", label="BVH model small")
plt.bar(
    r1 + 0.5,
    execution_times_long,
    width=0.5,
    color="green",
    label="AABB model small",
)
plt.bar(
    r1 + width1,
    execution_times_large,
    width=0.5,
    color="red",
    label="BVH model large",
)
plt.legend()
plt.xlabel("trial")
plt.ylabel("Execution Time")
plt.title("Execution Time bar graph")
plt.savefig("./bargraph.png")


def writeInfo(arr, f):
    f.write("Maximum execution time: %f \n" % (arr[len(arr) - 1]))
    f.write("Minimum execution time: %f \n" % (arr[0]))
    f.write("Mean execution time: %f \n" % (sum(arr) / len(arr)))

    if len(arr) % 2 == 0:
        f.write(
            "Median execution time: %f \n"
            % ((arr[len(arr) // 2] + arr[len(arr) // 2 + 1]) / 2)
        )
    else:
        f.write("Median execution time: %f \n" % (arr[len(arr) // 2]))


file = open("./result.txt", "w+")
execution_times_long.sort()

file.write("Execution time analysis for AABB structure small dataset \n")
writeInfo(execution_times_long, file)
execution_times_short.sort()
file.write("Execution time analysis for BVH structure small dataset \n")
writeInfo(execution_times_short, file)
execution_times_large.sort()
file.write("Execution time analysis for BVH structure large dataset \n")
writeInfo(execution_times_large, file)

file.close()
