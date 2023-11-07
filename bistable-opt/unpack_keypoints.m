function [rA1, rB1, rC1, rD1, rE1, rF1, rA2, rB2, rC2, rD2, rE2, rF2] = unpack_keypoints(keypoints)

    rA1 = keypoints(1:2,1);
    rB1 = keypoints(1:2,2);
    rC1 = keypoints(1:2,3);
    rD1 = keypoints(1:2,4);
    rE1 = keypoints(1:2,5);
    rF1 = keypoints(1:2,6);

    rA2 = keypoints(3:4,1);
    rB2 = keypoints(3:4,2);
    rC2 = keypoints(3:4,3);
    rD2 = keypoints(3:4,4);
    rE2 = keypoints(3:4,5);
    rF2 = keypoints(3:4,6);

end