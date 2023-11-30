import sys
import matplotlib.pyplot as plt
import math

# Get the array from command-line arguments
time_points = sys.argv[1:]

# Now you can use the execution_times array in your Python script
for i in range(len(time_points)):
    time_points[i] = float(time_points[i])

execution_times_short = []
execution_times_long = []
xaxis = []

for i in range(0, len(time_points), 2):
    exe = round((time_points[i + 1] - time_points[i]), 4)
    execution_times.append(exe)
    xaxis.append(i / 2 + 1)

fig = plt.figure(figsize=(10, 15))

# creating the bar plot
plt.bar(xaxis, execution_times, color="maroon", width=0.4)

plt.xlabel("trial")
plt.ylabel("Execution Time")
plt.title("Execution Time bar graph")
plt.savefig("./bargraph.png")

print(execution_times)
execution_times.sort()
print(execution_times)
