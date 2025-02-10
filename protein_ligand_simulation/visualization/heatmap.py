import seaborn as sns
import matplotlib.pyplot as plt

def generate_heatmap(data):
    sns.heatmap(data, annot=True, cmap="coolwarm")
    plt.show()
