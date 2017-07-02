from kinematics import *


def testsym():
    joint_angles = [-2.48, -0.76, 0.78, -0.68, -1.46, 3.63]

    print(calculate_FK(joint_angles))



if __name__ == "__main__":
    testsym()