from kinematics import *
from geometry_msgs.msg import Pose


def testk():
    pose = Pose()
    # pose.position.x = 2.1529
    # pose.position.z = 1.9465

    # pose.position.x = 1.8573
    # pose.position.y = 1.0887
    # pose.position.z = 1.9465


    # pose.orientation.x = 3.886783e-05
    # pose.orientation.y = -0.000143171
    # pose.orientation.z = 0.261995
    # pose.orientation.w = 0.965069

    # 0.53, 0, 0, -2.55, 0, 0
    # pose.orientation.x = 0.92362
    # pose.orientation.y = 0.25079
    # pose.orientation.z = -0.075807
    # pose.orientation.w = -0.27978

    # 0.53, 0.32, 0, -2.55, 0, 0
    # pose.position.x = 2.1008
    # pose.position.y = 1.2314
    # pose.position.z = 1.3283
    # pose.orientation.x = 0.924112
    # pose.orientation.y = 0.203753
    # pose.orientation.z = -0.219914
    # pose.orientation.w = -0.236937

    # 0.53, 0.32, -0.60, -2.55, 0, 0
    # pose.position.x = 2.1413
    # pose.position.y = 1.2551
    # pose.position.z = 2.3975
    # pose.orientation.x = 0.903253
    # pose.orientation.y = 0.288208
    # pose.orientation.z = 0.05705
    # pose.orientation.w = -0.312754

    # 0, 0, 0, 4.01, 0.99, 2.93
    # pose.position.x = 2.0153
    # pose.position.y = -0.1941
    # pose.position.z = 2.1101
    # pose.orientation.x = 0.2844
    # pose.orientation.y = -0.408428
    # pose.orientation.z = -0.245389
    # pose.orientation.w = 0.831981

    # 0, 0, 0, 0, 0, 0.5
    # pose.position.x = 2.1529
    # pose.position.y = 0
    # pose.position.z = 1.9465
    # pose.orientation.x = 0.247845
    # pose.orientation.y = -0.000143724
    # pose.orientation.z = 0
    # pose.orientation.w = 0.9688

    tests = [
        [
            [2.153, 0.000, 1.947],
            [0.000, 0.000, 0.000, 1.000],
            [0, 0, 0, 0, 0, 0]
        ],
        [
            [1.857, 1.089, 1.947],
            [0.000, -0.000, 0.262, 0.965],
            [0.53, 0, 0, 0, 0, 0]
        ],
        [
            [-0.693, -0.784, 1.805],
            [0.593, -0.587, 0.113, 0.540],
            [-2.48, -0.76, 0.78, -0.68, -1.46, 3.63]
        ]
    ]


    for test in tests:
        position = test[0]
        orientation = test[1]
        joint_angles = test[2]

        pose.position.x = position[0]
        pose.position.y = position[1]
        pose.position.z = position[2]
        pose.orientation.x = orientation[0]
        pose.orientation.y = orientation[1]
        pose.orientation.z = orientation[2]
        pose.orientation.w = orientation[3]

        print ("----")
        check_IK(pose, calculate_IK(pose))
        check_FK(pose, calculate_FK(joint_angles))

if __name__ == "__main__":
    testk()