## Project: Kinematics Pick & Place
### Craig Robinson

---

[//]: # (Image References)

[theta3]: ./misc_images/IMG_6701.JPG
[theta2]: ./misc_images/IMG_6702.JPG
[dhparams]: ./misc_images/IMG_6703.JPG
[complete_run]: ./misc_images/complete_run.png

## [Rubric](https://review.udacity.com/#!/rubrics/972/view) Points
### Here I will consider the rubric points individually and describe how I addressed each point in my implementation.


### Kinematic Analysis
#### 1. 
I used the following diagram to determine DH parameters:

![DH Parameter Diagram][dhparams]

Here is my DH parameter table:

    s = {
        alpha0:       0, a0:      0, d1:  0.75,                   # 0 -> 1
        alpha1: -pi / 2, a1:   0.35, d2:     0,  q2: q2 - pi / 2, # 1 -> 2
        alpha2:       0, a2:   1.25, d3:     0,                   # 2 -> 3
        alpha3: -pi / 2, a3: -0.054, d4:  1.50,                   # 3 -> 4
        alpha4:  pi / 2, a4:      0, d5:     0,                   # 4 -> 5
        alpha5: -pi / 2, a5:      0, d6:     0,                   # 5 -> 6
        alpha6:       0, a6:      0, d7: 0.303, q7: 0
    }

This is how each one of the parameters was defined from the URDF file:

    alpha0 = 0       # no twist Z0 -> Z1
    alpha1 = -pi / 2 # -90 degree twist Z1 -> Z2
    alpha2 = 0       # parallel Z2 -> Z3
    alpha3 = -pi / 2 # -90 degree twist Z3 -> Z4
    alpha4 = pi / 2  # 90 degree twist Z4 -> Z5
    alpha5 = -p1 / 2 # -90 degree twist Z5 -> Z6
    alpha6 = 0       # parallel Z6 -> gripper
    a0 = 0           # base
    a1 = 0.35        # Z1 -> Z2 along X1 axis: joint_2.x
    a2 = 1.25        # Z2 -> Z3 along X2 axis: joint_3.z
    a3 = -0.054      # Z3 -> Z4 along X3: joint_4.z
    a4 = 0           # O4, O5 are coincident
    a5 = 0           # O5, O6 are coincident
    a6 = 0           # Z6, ZG are colinear
    d1 = 0.75        # X0 -> X1 along Z1 axis: joint_1.z + joint_2.z    
    d2 = 0 #         # O1, O2 are coincident
    d3 = 0           # X2, X3 are colinear
    d4 = 1.50        # X3 -> X4 along Z4: joint_4.x + joint_5.x
    d6 = 0           # X5, X6 are colinear
    d7 = 0.303       # distance from 06 to O_gripper: joint_6.x + gripper.x

#### 2. Matrices

Individual transform matrices about each joint:

    ('T0_1 = ', Matrix([
    [cos(q1), -sin(q1), 0,    0],
    [sin(q1),  cos(q1), 0,    0],
    [      0,        0, 1, 0.75],
    [      0,        0, 0,    1]]))
    
    ('T1_2 = ', Matrix([
    [ cos(q2), -sin(q2), 0, 0.35],
    [       0,        0, 1,    0],
    [-sin(q2), -cos(q2), 0,    0],
    [       0,        0, 0,    1]]))
    
    ('T2_3 = ', Matrix([
    [cos(q3), -sin(q3), 0, 1.25],
    [sin(q3),  cos(q3), 0,    0],
    [      0,        0, 1,    0],
    [      0,        0, 0,    1]]))
    
    ('T3_4 = ', Matrix([
    [ cos(q4), -sin(q4), 0, -0.054],
    [       0,        0, 1,    1.5],
    [-sin(q4), -cos(q4), 0,      0],
    [       0,        0, 0,      1]]))
    
    ('T4_5 = ', Matrix([
    [cos(q5), -sin(q5),  0, 0],
    [      0,        0, -1, 0],
    [sin(q5),  cos(q5),  0, 0],
    [      0,        0,  0, 1]]))
    
    ('T5_6 = ', Matrix([
    [ cos(q6), -sin(q6), 0, 0],
    [       0,        0, 1, 0],
    [-sin(q6), -cos(q6), 0, 0],
    [       0,        0, 0, 1]]))
    
Each matrix is formed by substituting the DH parameters into the DH transform matrix. See calculate_FK() in kinematics.py.

Generalized homogeneous transform between base_link and gripper_link using only end-effector(gripper) pose:

    ⎡sin(pitch)⋅cos(roll)⋅cos(yaw) + sin(roll)⋅sin(yaw)  -sin(pitch)⋅sin(roll)⋅cos(yaw) + sin(yaw)⋅cos(roll)  cos(pitch)⋅cos(yaw)  px⎤
    ⎢                                                                                                                                ⎥
    ⎢sin(pitch)⋅sin(yaw)⋅cos(roll) - sin(roll)⋅cos(yaw)  -sin(pitch)⋅sin(roll)⋅sin(yaw) - cos(roll)⋅cos(yaw)  sin(yaw)⋅cos(pitch)  py⎥
    ⎢                                                                                                                                ⎥
    ⎢               cos(pitch)⋅cos(roll)                                -sin(roll)⋅cos(pitch)                     -sin(pitch)      pz⎥
    ⎢                                                                                                                                ⎥
    ⎣                        0                                                    0                                    0           1 ⎦

This matrix was created by the composition of the roll, pitch and yaw values (converted from a quaternion expressing the 
pose of the gripper), followed by a correction matrix to convert from URDF to the base coordinate frame, and finally, 
by adding in the translation from the gripper pose. See transform_base_gripper() in kinematics.py.

#### 3. Decoupled Inverse Kinematics

##### a: Inverse Position

###### i. Computing wrist center

To compute the wrist center I first computed R_corr to convert from URDF to DH space. It is a rotation about the Z axis 
by 90 degrees followed by a rotation about the Y axis by -90 degrees. I then compute the rotation matrix for the wrist
by using the roll, pitch, and yaw values of the end effector pose (converted from the quarternion). It is an extrinsic 
rotation sequence `Rrpy = Rz(yaw) * Ry(pitch) * Rx(roll) * R_corr`. Finally, I find the wrist center by translation along
the Z axis of the gripper the distance from the gripper to the wrist center.

    # correction from URDF to DH (zy intrinsic)
    R_corr = rot_z(pi) * rot_y(-pi / 2)

    # compute rotation matrix for gripper (zyx intrinsic / xyz extrinsic)
    Rrpy = rot_z(yaw) * rot_y(pitch) * rot_x(roll) * R_corr.transpose()

    # wrist center = p - (d6 + leef) * Znormal
    wx = px - (d6 + d7) * Rrpy[0, 2]
    wy = py - (d6 + d7) * Rrpy[1, 2]
    wz = pz - (d6 + d7) * Rrpy[2, 2]

    
###### i. Computing theta3

I determined the equations for joint3 by isolating the triangle formed by links 2 and 3 in the X-Z plane and using the 
cosine rule to find the joint angle. 

The distal end of link3, O4, is at wrist center (as computed above). The proximal end of link2 is at O2, which is at a1, 
d1 in the X-Z plane.

The generalized location of the wrist center, projected onto the X-Y plane can be computed as `wxy = sqrt(wx * wx + wy * wy)`. We
can use this value and the location of the wrist center in the z direction, wz, to compute the length of the vector from
O2 to O4 (WC), as `c = sqrt((wxy - a1) * (wxy - a1) + (wz - d1) * (wz - d1))`.

The side of the triangle represented by link 2 has a length of a2.

The side of the triangle represented by link 3 has a length of `b = sqrt(a3 * a3 + d4 * d4)`, since there is a small offset
between the z location of joint3 and joint4.

Given the lengths of the sides of the triangle, we can apply the cosine rule to find the cosine and sine of the angle between
links 2 and 3 as `cos_theta3b = (a**2 + b**2 - c**2) / (2 * a * b)` and `sin_theta3b = sqrt(1 - cos_theta3b**2)`. Finally,
we use the atan2 function to find the angle between links 2 and 3 and then subtract the fixed offset of the joint,
 `theta3 = atan2(sin_theta3, -cos_theta3) - atan2(d4, a3)`

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
    
![Computing Theta3][theta3]

###### ii. Computing theta2

Theta2 is computed by extending the line formed by link2 beyond joint3 to form a right triangle with b as the hypotenuse.
The supplementary angle to theta3b (theta2c = 180 - theta3b) is one of the angles of this triangle. We can use these values
to determine the angle theta2b. Theta2a is the angle between the X axis and the line from 02 to the wrist center. 
Then theta2 = 90 - (theta2a + theta2b).

    # theta2
    theta2a = atan2(wz - d1, wxy - a1)
    theta2c = pi - theta3b
    cos_theta2c = cos(theta2c)
    sin_theta2c = sin(theta2c)
    theta2b = atan2(b * sin_theta2c, a2 + -b * cos_theta2c)
    theta2 = pi / 2 - (theta2a + theta2b)
    
![Computing Theta2][theta2]
     

###### iii. Computing theta1

Theta1 is easily computed as the atan of of the distance to the wrist center in the Y and X axes.

    theta1 = atan2(wy, wx)
    

##### b: Inverse Orientation

The inverse orientation problem was solved by computing `R3_6 = R0_3.transpose() * Rrpy`. The trick to this was to first
compute R3_6 symbolically and analyze the resulting rotation matrix to determine equations for the joint angles. I used
a method similar to the lesson on Euler angles from a rotation matrix to isolate the sine & cosine values for the joints 
and use atan2 to compute the angles. One tricky part was handling the singularity case. I check to see if the cos(q5) is 
close to 0, and if so, set q5 to 90, q5 to 0 and compute q6 from values in the matrix.

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
 

### Project Implementation

#### 1. Code

My code successfully completes the pick and place problem. In this run it successfully moved 9 out of 10 pegs. The one
that it failed on looked correct from an IK standpoint, but for some reason, it failed to grasp the peg.

The pose error rates are on the order of 1.0e-16 meters/radians.

![Complete Run][complete_run]

I implemented most of the code in a file called kinematics.py. I call this from IK_server.py for the simulation, but
I also created a couple of other top-level scripts to call the functions in kinematics.py for testing.

I used the recommended techniques from the lesson to solve the IK problem. I used SymPy for the forward kinematics,
but found it difficult to use (mostly due to my unfamiliarity with it) so I solved the IK problem without it. I think
I understand now how I could use SymPy in the IK problem and it might be interesting to go back and do that now that 
I solved the problem without it.

One of the things I noticed with my implementation is that the wrist spins more than it needs to. An improvement would
be to go back and figure out how to minimize wrist movement, perhaps by keeping track of the current angle and computing
the minimum turn direction needed to achieve the new orientation.


