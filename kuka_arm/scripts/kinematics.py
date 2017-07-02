
import tf
from mpmath import *
from sympy import *

def within_eps(v, t, err):
    return abs(v) > (t - err) and abs(v) < (t + err)

def within(v, t):
    eps = 1.0e-8
    return within_eps(v, t, eps)

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


def rot2euler(R):
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

def calculate_FK(joint_angles):
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

    if len(joint_angles) >= 6:
        q1 = q1.subs({q1: joint_angles[0]})
        q2 = q2.subs(s).subs({q2: joint_angles[1]})
        q3 = q3.subs({q3: joint_angles[2]})
        q4 = q4.subs({q4: joint_angles[3]})
        q5 = q5.subs({q5: joint_angles[4]})
        q6 = q6.subs({q6: joint_angles[5]})
    else:
        print ("Solving symbolically")

    T0_1 = dht(alpha0, a0, d1, q1)
    T0_1 = T0_1.subs(s)
    T1_2 = dht(alpha1, a1, d2, q2)
    T1_2 = T1_2.subs(s)
    T2_3 = dht(alpha2, a2, d3, q3)
    T2_3 = T2_3.subs(s)
    T3_4 = dht(alpha3, a3, d4, q4)
    T3_4 = T3_4.subs(s)
    T4_5 = dht(alpha4, a4, d5, q5)
    T4_5 = T4_5.subs(s)
    T5_6 = dht(alpha5, a5, d6, q6)
    T5_6 = T5_6.subs(s)
    T6_G = dht(alpha6, a6, d7, q7)
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

    T_total = simplify(T0_G * R_corr)

    return T_total

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


    # print ("r23", r23)
    if within_eps(r23, 0, 0.001): # q5 ~= 0 => singularity
        # print ("SINGULARITY")
        q5 = pi

        q4 = Float(0)

        # r12 = -sin(q4) * cos(q6) - sin(q6) * cos(q4) * cos(q5)
        # = -sin(q4) * 1 - 0 * cos(q4) * cos(q5)
        # = -sin(q4)

        # r32 = sin(q4) * sin(q6) * cos(q5) - cos(q4) * cos(q6)
        # = sin(q4) * 0 * cos(q5) - cos(q4) * 1
        # = -cos(q4)

        q6 = atan2(-r22, r21)

    else:
        q5 = atan2(sqrt((r13 * r13 + r33 * r33)), r23)

        # r22 = -sin(q5) * sin(q6)
        # r21 = sin(q5) * cos(q6)
        # -r22/r21 = sin(q6)/cos(q6)
        q6 = atan2(-r22, r21)

        # r33 = sin(q4) * sin(q5)
        # r13 = -sin(q5) * cos(q4)
        # r33/-r13 = sin(q4)/cos(q4)

        q4 = atan2(r33, -r13)

    if within(q4, pi):
        q4 = Float(0)

    if within(q5, pi):
        q5 = Float(0)

    if within(q6, pi):
        q6 = Float(0)

    return q4, q5, q6


def transform_base_gripper():

    px = symbols('px')
    py = symbols('py')
    pz = symbols('pz')
    roll = symbols('roll')
    pitch = symbols('pitch')
    yaw = symbols('yaw')

    R_corr = rot_z(pi) * rot_y(-pi / 2)
    Rrpy = rot_z(yaw) * rot_y(pitch) * rot_x(roll) * R_corr.transpose()

    T = Rrpy.row_join(Matrix([[px], [py], [pz]]))
    T = T.col_join(Matrix([[0, 0, 0, 1]]))

    pprint(T)


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

    # correction from URDF to DH (zy intrinsic)
    R_corr = rot_z(pi) * rot_y(-pi / 2)

    # compute rotation matrix for gripper (zyx intrinsic / xyz extrinsic)
    Rrpy = rot_z(yaw) * rot_y(pitch) * rot_x(roll) * R_corr.transpose()

    # wrist center = p - (d6 + leef) * Znormal
    wx = px - (d6 + d7) * Rrpy[0, 2]
    wy = py - (d6 + d7) * Rrpy[1, 2]
    wz = pz - (d6 + d7) * Rrpy[2, 2]

    # print("w", wx, wy, wz)

    # temp vars for computing theta3 & theta2
    a = a2
    asqr = a**2
    bsqr = a3 * a3 + d4 * d4
    b = sqrt(bsqr)
    wxy = sqrt(wx**2 + wy**2)
    csqr = (wxy - a1)**2 + (wz - d1)**2
    c = sqrt(csqr)

    # theta 3
    theta3a = atan2(d4, a3)
    cos_theta3b = (asqr + bsqr - csqr) / (2 * a * b) # cosine rule
    sin_theta3b = sqrt(1 - cos_theta3b**2)           # sin^2 + cos^2 = 1
    theta3b = atan2(sin_theta3b, -cos_theta3b)
    theta3 = theta3b - theta3a

    # theta2
    theta2a = atan2(wz - d1, wxy - a1)
    theta2c = pi - theta3b
    cos_theta2c = cos(theta2c)
    sin_theta2c = sin(theta2c)
    theta2b = atan2(b * sin_theta2c, a2 + -b * cos_theta2c)
    theta2 = pi / 2 - (theta2a + theta2b)

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
    (yawp, pitchp, rollp) = rot2euler(T_check)

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

    print "trans_error: %s" % error_bar(avg_trans_error)
    print "  rot_error: %s" % error_bar(avg_rot_error)

    if avg_trans_error > eps or avg_rot_error > eps:
        print("******************")
        print("* ERROR TOO HIGH *")
        print("******************")

def error_bar(error):

    error_bar = ""
    while error > 1.0e-16:
        error_bar += "*"
        error /= 10

    return error_bar