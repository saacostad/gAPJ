import numpy as np
import matplotlib.pyplot as plt
import subprocess
import pandas as pd

N = 500
error_counter = 0


df = pd.DataFrame(columns=["x_before", "y_before", "x_after", "y_after"])



for i in range(N):

    print(f"\n========== Run {i + 1}/{N} ==========\n")

    try:
        result = subprocess.run(
            ["bash", "commands"],
            capture_output=True,
            text=True,
            check=True
        )

        last_line = result.stdout.splitlines()[-1]

        print(f"Last output: {last_line}")

        df.loc[len(df)] = last_line

    except subprocess.CalledProcessError as e:
        print(f"\nPipeline failed on iteration {i + 1}")
        print(f"Return code: {e.returncode}")
        error_counter += 1
        continue

    else:
        print(f"\nIteration {i + 1} completed successfully!")
    

print("="*20)
print("FINISHED HISTOGRAM")
print("# of simulations: ", N - error_counter)
df.to_csv("LHC.csv", index=False)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# X histogram
ax[0].hist(df["x_before"], bins = 30, alpha=0.6, label="Before")
ax[0].hist(df["x_after"], bins = 30, alpha=0.6, label="After")
ax[0].set_xlabel("X")
ax[0].set_ylabel("Frequency")
ax[0].set_title(r"X $\beta$-beating distribution")
ax[0].legend()

# Y histogram
ax[1].hist(df["y_before"], bins = 30, alpha=0.6, label="Before")
ax[1].hist(df["y_after"], bins = 30, alpha=0.6, label="After")
ax[1].set_xlabel("Y")
ax[1].set_ylabel("Frequency")
ax[1].set_title(r"Y $\beta$-beating distribution")
ax[1].legend()

plt.tight_layout()
plt.show()
