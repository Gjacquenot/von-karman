import h5py

def save_step(folder_path, xx, yy, p, u, v, t):
    with h5py.File(folder_path+f"U_{round(t,4)}", "w") as f:
        f.create_dataset("X", data=xx)
        f.create_dataset("Y", data=yy)
        f.create_dataset("p", data=p)
        f.create_dataset("u", data=u)
        f.create_dataset("v", data=v)
        f.create_dataset("t", data=t)
    return None