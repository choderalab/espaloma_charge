import pandas as pd
import numpy as np

def bootstrap_mean(x, n_samples=1000, ci=0.95):
    results = []
    for _ in range(n_samples):
        idxs = np.random.randint(x.shape[0], size=(x.shape[0], ))
        results.append(x[idxs].mean())

    low = np.percentile(results, 100.0 * 0.5 * (1 - ci))
    high = np.percentile(results, (1 - ((1 - ci) * 0.5)) * 100.0)
    original = x.mean()
    return original, low, high
        




def run(in_path):
    df = pd.read_csv(
        in_path, 
        names=["smiles", "q_amber", "q_esp", "q_oe", "t_amber", "t_esp", "t_oe"],
    )

    for name in ["q_amber", "q_esp", "q_oe"]:
        df[name] = df.apply(lambda x: np.array(eval(x[name])), axis=1)

    df["q_rmse_esp"] = ((df["q_esp"] - df["q_oe"]) ** 2).apply(np.mean) ** 0.5
    df["q_rmse_amber"] = ((df["q_amber"] - df["q_oe"]) ** 2).apply(np.mean) ** 0.5
    df["n_atoms"] = df.apply(lambda x: len(x["q_oe"]), axis=1)

    n_mol = "$ %s $" % len(df)
    n_atoms = "$ %.2f $" % df["n_atoms"].mean() 

    q_rmse_esp = "$%.4f_{%.4f}^{%.4f}$" % bootstrap_mean(df["q_rmse_esp"].values.flatten())
    q_rmse_amber = "$%.4f_{%.4f}^{%.4f}$" % bootstrap_mean(df["q_rmse_amber"].values.flatten())

    t_esp = "$ %.2f $" % df["t_esp"].mean()
    t_amber = "$ %.2f $" % df["t_amber"].mean()
    t_oe = "$ %.2f $" % df["t_oe"].mean()

    print(" & ".join([n_mol, n_atoms, q_rmse_esp, q_rmse_amber, t_esp, t_amber, t_oe]))

if __name__ == "__main__":
    import sys
    run(sys.argv[1])

    
