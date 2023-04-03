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

with open(os.path.join(root, "Tests.m"), "w") as f:
    print("Tests");
    path = os.path.join(root, "Tests")
    f.write("AttachSpec(\"spec\");\n")
    for test in [test for test in os.listdir(path) if test[-1] == "m"]:
        f.write("load \"")
        f.write(os.path.join(path, test));
        f.write("\";\n")
        f.write("print(\"PASS\");\n");
    f.write("exit;");
