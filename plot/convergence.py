import pandas as pd
import sys
import matplotlib.pyplot as plt

#df = pd.read_csv("convergence.csv", delimiter=" ")
df_gs  = pd.read_csv("gs.csv", delimiter=" ")
df_sor = pd.read_csv("sor.csv", delimiter=" ")

fig, axes = plt.subplots()

#df.plot(ax=axes, x="time", y="error", label="Current")
df_sor.plot(ax=axes, x="time", y="error", label="SOR", color="red")
df_gs.plot( ax=axes, x="time", y="error", label="GS", color="blue")

plt.grid(True)
plt.title("Heat 2D steady state Gauss-Seidel vs Succesive-Over-Relaxation")
plt.ylabel("Absolute error (K)")
plt.xlabel("Time (s)")
plt.yscale("log")
plt.savefig("err.png")

