import os
import matplotlib.pyplot as plt
import pandas as pd

from src.logging import lprint, LoggingLevels as ll


def save_distribution_plot(series: pd.Series, folder: str, target: str, bins: int = 30):
    os.makedirs(folder, exist_ok=True)
    
    plt.figure(figsize=(8, 5))
    plt.hist(series, bins=bins, density=True, alpha=0.6, edgecolor="black")
    series.plot(kind="kde", color="red", label="KDE")
    plt.title("Distribution")
    plt.xlabel("Valore")
    plt.ylabel("Densità")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    filepath = os.path.join(folder, f"distribution_{target}.png")
    plt.savefig(filepath, dpi=300, bbox_inches="tight")
    plt.close()  
    lprint(ll.REPORT, f"Distribution plot saved in {filepath}")
