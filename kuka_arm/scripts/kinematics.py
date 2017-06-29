
from geometry_msgs.msg import Pose
import tf
from mpmath import *
from sympy import *
from sympy.printing.str import StrPrinter

def dtr(d):
    return (d * pi / 180.0).evalf()

def rtd(r):
    return (r * 180.0 / pi).evalf()

### Define functions for Rotation Matrices about x, y, and z given specific angle.
def rot_x(q):
    R_x = Matrix([[1, 0, 0],
                  [0, cos(q), -sin(q)],
                  [0, sin(q), cos(q)]])

    return R_x


def rot_y(q):
    R_y = Matrix([[cos(q), 0, sin(q)],
                  [0, 1, 0],
                  [-sin(q), 0, cos(q)]])

    return R_y


def rot_z(q):
    R_z = Matrix([[cos(q), -sin(q), 0],
                  [sin(q), cos(q), 0],
                  [0, 0, 1]])

    return R_z


def dht(alpha, a, d, q):
    dh = Matrix([[cos(q), -sin(q), 0, a],
                 [sin(q) * cos(alpha), cos(q) * cos(alpha), -sin(alpha), -sin(alpha) * d],
                 [sin(q) * sin(alpha), cos(q) * sin(alpha), cos(alpha), cos(alpha) * d],
                 [0, 0, 0, 1]])

    return dh


def rot2euler_orig(R):
    r12 = R[0, 1]
    r13 = R[0, 2]
    r31 = R[2, 0]
    r11 = R[0, 0]
    r21 = R[1, 0]
    r32 = R[2, 1]
    r33 = R[2, 2]

    cos_beta = sqrt(r11 * r11 + r21 * r21)
    beta = atan2(-r31, cos_beta)
    gamma = atan2(r32, r33)
    alpha = atan2(r21, r11)

    return alpha, beta, gamma

def symbolic_FK():
    q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
    d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
    a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
    alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

    # Modified DH params
    s = {
        alpha0: 0, a0: 0, d1: 0.75,  # 0 -> 1
        alpha1: -pi / 2, a1: 0.35, d2: 0,  q2: q2 - pi / 2, # 1 -> 2
        alpha2: 0, a2: 1.25, d3: 0,  # 2 -> 3
        alpha3: -pi / 2, a3: -0.054, d4: 1.50,  # 3 -> 4
        alpha4: pi / 2, a4: 0, d5: 0,  # 4 -> 5
        alpha5: -pi / 2, a5: 0, d6: 0,  # 5 -> 6
        alpha6: 0, a6: 0, d7: 0.303, q7: 0
    }

    T0_1 = dht(alpha0, a0, d1, q1)
    T0_1.subs(s)
    T1_2 = dht(alpha1, a1, d2, q2)
    T1_2.subs(s)
    T2_3 = dht(alpha2, a2, d3, q3)
    T2_3.subs(s)
    T3_4 = dht(alpha3, a3, d4, q4)
    T3_4.subs(s)
    T4_5 = dht(alpha4, a4, d5, q5)
    T4_5.subs(s)
    T5_6 = dht(alpha5, a5, d6, q6)
    T5_6.subs(s)
    T6_G = dht(alpha6, a6, d7, q7)
    T6_G.subs(s)

    T0_3 = simplify(T0_1 * T1_2 * T2_3)
    T0_3 = T0_3.subs(s)
    T3_6 = simplify(T3_4 * T4_5 * T5_6)
    T3_6 = T3_6.subs(s)

    R_z = Matrix([[cos(pi), -sin(pi), 0, 0],
                  [sin(pi), cos(pi), 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])

    R_y = Matrix([[cos(-pi / 2), 0, sin(-pi / 2), 0],
                  [0, 1, 0, 0],
                  [-sin(-pi / 2), 0, cos(-pi / 2), 0],
                  [0, 0, 0, 1]])

    R_corr = simplify(R_z * R_y)

    T0_G = simplify(T0_3 * T3_6 * T6_G * R_corr)
    T0_G.subs(s)

    # print ("T3_6", T3_6)
    # print ("T0_G", T0_G)

    return T0_G

def joint_angles_from_3_6(R3_6):

    # derive equations for the last 3 joint angles from R3_6:
    # [-sin(q4) * sin(q6) + cos(q4) * cos(q5) * cos(q6), -sin(q4) * cos(q6) - sin(q6) * cos(q4) * cos(q5), -sin(q5) * cos(q4)],
    # [                               sin(q5) * cos(q6),                               -sin(q5) * sin(q6),            cos(q5)],
    # [-sin(q4) * cos(q5) * cos(q6) - sin(q6) * cos(q4),  sin(q4) * sin(q6) * cos(q5) - cos(q4) * cos(q6),  sin(q4) * sin(q5)],

    r12 = R3_6[0, 1] # -sin(q4) * cos(q6) - sin(q6) * cos(q4) * cos(q5)
    r13 = R3_6[0, 2] # -sin(q5) * cos(q4)
    r21 = R3_6[1, 0] # sin(q5) * cos(q6)
    r22 = R3_6[1, 1] # -sin(q5) * sin(q6)
    r23 = R3_6[1, 2] # cos(q5)
    r32 = R3_6[2, 1] # sin(q4) * sin(q6) * cos(q5) - cos(q4) * cos(q6)
    r33 = R3_6[2, 2] # sin(q4) * sin(q5)

    # r13 * r13 + r33 * r33
    # = (-sin(q5) * cos(q4)) * (-sin(q5) * cos(q4)) + (sin(q4) * sin(q5)) + (sin(q4) * sin(q5))
    # = sin(q5)^2 * cos(q4)^2 + sin(q4)^2 * sin(q5)2
    # = sin(q5)^2 * (cos(q4)^2 + sin(q4)^2)
    # = sin(q5)^2
    q5 = atan2(sqrt((r13 * r13 + r33 * r33)), r23)

    if abs(q5) < 0.0e-5: # q5 ~= 0 => singularity

        q6 = Float(0)

        # r12 = -sin(q4) * cos(q6) - sin(q6) * cos(q4) * cos(q5)
        # = -sin(q4) * 1 - 0 * cos(q4) * cos(q5)
        # = -sin(q4)

        # r32 = sin(q4) * sin(q6) * cos(q5) - cos(q4) * cos(q6)
        # = sin(q4) * 0 * cos(q5) - cos(q4) * 1
        # = -cos(q4)

        q4 = atan2(r12, r32)

    else:


        # r22 = -sin(q5) * sin(q6)
        # r21 = sin(q5) * cos(q6)
        # -r22/r21 = sin(q6)/cos(q6)
        q6 = atan2(-r22, r21)

        # r33 = sin(q4) * sin(q5)
        # r13 = -sin(q5) * cos(q4)
        # r33/-r13 = sin(q4)/cos(q4)

        q4 = atan2(r33, -r13)

        return q4, q5, q6



def calculate_FK(joint_angles):
    q1 = joint_angles[0]
    q2 = joint_angles[1]
    q3 = joint_angles[2]
    q4 = joint_angles[3]
    q5 = joint_angles[4]
    q6 = joint_angles[5]

    d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
    a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
    alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

    # Modified DH params
    s = {
        alpha0: 0, a0: 0, d1: 0.75,  # 0 -> 1
        alpha1: -pi / 2, a1: 0.35, d2: 0,  # 1 -> 2
        alpha2: 0, a2: 1.25, d3: 0,  # 2 -> 3
        alpha3: -pi / 2, a3: -0.054, d4: 1.50,  # 3 -> 4
        alpha4: pi / 2, a4: 0, d5: 0,  # 4 -> 5
        alpha5: -pi / 2, a5: 0, d6: 0,  # 5 -> 6
        alpha6: 0, a6: 0, d7: 0.303
    }

    T0_1 = dht(alpha0, a0, d1, q1)
    T0_1 = T0_1.subs(s)
    T1_2 = dht(alpha1, a1, d2, q2 - pi / 2)
    T1_2 = T1_2.subs(s)
    T2_3 = dht(alpha2, a2, d3, q3)
    T2_3 = T2_3.subs(s)
    T3_4 = dht(alpha3, a3, d4, q4)
    T3_4 = T3_4.subs(s)
    T4_5 = dht(alpha4, a4, d5, q5)
    T4_5 = T4_5.subs(s)
    T5_6 = dht(alpha5, a5, d6, q6)
    T5_6 = T5_6.subs(s)
    T6_G = dht(alpha6, a6, d7, 0)
    T6_G = T6_G.subs(s)

    T0_G = simplify(T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_G)

    R_z = Matrix([[cos(pi), -sin(pi), 0, 0],
                  [sin(pi), cos(pi), 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])

    R_y = Matrix([[cos(-pi / 2), 0, sin(-pi / 2), 0],
                  [0, 1, 0, 0],
                  [-sin(-pi / 2), 0, cos(-pi / 2), 0],
                  [0, 0, 0, 1]])

    R_corr = simplify(R_z * R_y)
    T0_G = T0_G * R_corr

    T_total = simplify(T0_G)

    return T_total


def calculate_IK(pose):

    # dh params
    alpha0 = 0
    alpha1 = -pi / 2
    alpha2 = 0
    a0 = 0
    a1 = 0.35 # joint_2.x
    d1 = 0.75 # joint_1.z + joint_2.z
    a2 = 1.25 # joint_3.z
    d2 = 0
    a3 = -0.054 # joint_4.z
    d3 = 0
    d4 = 1.50 # joint_4.x + joint_5.x
    d6 = 0
    d7 = 0.303 # joint_6.x + gripper.x

    # Extract end-effector position and orientation from request
    # px,py,pz = end-effector position
    # roll, pitch, yaw = end-effector orientation
    px = pose.position.x
    py = pose.position.y
    pz = pose.position.z

    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
        [pose.orientation.x, pose.orientation.y,
         pose.orientation.z, pose.orientation.w])

    # print ("roll", roll, "pitch", pitch, "yaw", yaw)
    # Calculate joint angles using Geometric IK method

    # compute rotation matrix for gripper
    R_corr = rot_z(pi) * rot_y(-pi / 2)
    Rrpy = rot_z(yaw) * rot_y(pitch) * rot_x(roll) * R_corr

    # wrist center = p - (d6 + leef) * Znormal
    wx = px - (d6 + d7) * Rrpy[0, 2]
    wy = py - (d6 + d7) * Rrpy[1, 2]
    wz = pz - (d6 + d7) * Rrpy[2, 2]

    # print("w", wx, wy, wz)

    # temp vars for computing theta3 & theta2
    a = a2
    b = sqrt(a3 * a3 + d4 * d4)
    wxy = sqrt(wx * wx + wy * wy)
    c2 = (wxy - a1) * (wxy - a1) + (wz - d1) * (wz - d1)

    # apply cosine rule to find theta 3
    cos_theta3 = (a * a + b * b - c2) / (2 * a * b)
    sin_theta3 = sqrt(1 - cos_theta3 * cos_theta3)
    theta3 = atan2(sin_theta3, -cos_theta3) - atan2(d4, a3)

    theta2 = pi / 2 - (atan2(wz - d1, wxy - a1) + atan2(b * sin_theta3, a2 + -b * cos_theta3))

    theta1 = atan2(wy, wx)

    # compute theta4-6
    T0_1 = dht(alpha0, a0, d1, theta1)
    T1_2 = dht(alpha1, a1, d2, theta2 - pi / 2)
    T2_3 = dht(alpha2, a2, d3, theta3)
    T0_3 = T0_1 * T1_2 * T2_3

    R0_3 = T0_3
    R0_3.row_del(3)
    R0_3.col_del(3)
    R3_6 = R0_3.transpose() * Rrpy
    (theta4, theta5, theta6) = joint_angles_from_3_6(R3_6)

    return [theta1.evalf(), theta2.evalf(), theta3.evalf(), theta4.evalf(), theta5.evalf(), theta6.evalf()]

def check_FK(pose, T):
    px = pose.position.x
    py = pose.position.y
    pz = pose.position.z
    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
        [pose.orientation.x, pose.orientation.y,
         pose.orientation.z, pose.orientation.w])

    Rp = rot_z(yaw) * rot_y(pitch) * rot_x(roll)
    Tp = Rp.row_join(Matrix([px, py, pz])).col_join(Matrix([[0, 0, 0, 1]]))

    # init_printing()
    # print ("\nTinput")
    # pprint (Tp)
    # print ("\nT_computed")
    # pprint (T)
    # print ("\nTerror")
    # pprint ((Tp - T))



def check_IK(pose, joint_angles):
    px = pose.position.x
    py = pose.position.y
    pz = pose.position.z

    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
        [pose.orientation.x, pose.orientation.y,
         pose.orientation.z, pose.orientation.w])

    T_check = calculate_FK(joint_angles)
    (yawp, pitchp, rollp) = rot2euler_orig(T_check)

    yawp = yawp.evalf()
    pitchp = pitchp.evalf()
    rollp = rollp.evalf()

    cpx = (T_check[0, 3]).evalf()
    cpy = (T_check[1, 3]).evalf()
    cpz = (T_check[2, 3]).evalf()
    Epx = px - cpx
    Epy = py - cpy
    Epz = pz - cpz
    Eroll = (roll - rollp).evalf()
    Epitch = (pitch - pitchp).evalf()
    Eyaw = (yaw - yawp).evalf()

    print("joint_angles: ", joint_angles)
    print("Input pose:    ", "px", px, "py", py, "pz", pz, "roll", roll, "pitch", pitch, "yaw", yaw)
    print("Computed pose: ", "px", cpx, "py", cpy, cpz, "roll", rollp, "pitch", pitchp, "yaw", yawp)
    print("Error          ", "px", Epx, "py", Epy, "pz", Epz, "roll", Eroll, "pitch", Epitch, "yaw", Eyaw)

    eps = 1.0e-10
    avg_trans_error = abs(Epx) + abs(Epy) + abs(Epz) / 3
    avg_rot_error = abs(Eroll) + abs(Epitch) + abs(Eyaw) / 3

    if avg_trans_error > eps or avg_rot_error > eps:
        print("******************")
        print("* ERROR TOO HIGH *")
        print("******************")
