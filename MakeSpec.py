import sys, os
root = "."

with open(os.path.join(root, "spec"), "w") as spec:
    spec.write("{\n")
    for D in [D for D in os.listdir(root) if os.path.isdir(D) and D[0] != "." and D != "Tests"]:
        print(D)
        #Generate spec file in subdirectory
        path = os.path.join(root, D)
        with open(os.path.join(path, "spec"), "w") as f:
            spec.write("    ")
            spec.write("+")
            spec.write(D)
            spec.write("/spec\n")
            f.write("{\n")
            for src in [src for src in os.listdir(path) if src[-1] == "m"]:
                f.write("    ")
                f.write(src)
                f.write("\n")
            f.write("}\n")
    spec.write("}\n")
