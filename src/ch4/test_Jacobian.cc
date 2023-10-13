//4.42 旋转残差对Ri的雅可比矩阵
    const SO3 R1 = p1->estimate().so3();
    const SO3 R1T = R1.inverse();
    const SO3 R2 = p2->estimate().so3();


    //获取零偏更新后的预积分测量值
    const SO3 dR = preint_->GetDeltaRotation(bg);
    const SO3 eR = SO3(dR).inverse() * R1T * R2;

    //自己修改的部分
    //思路是在算eR时，在R1基础上加右扰动
    //微小扰动量dd_phi
    const vector3d dd_phi;
    const SO3 eR_dd_phi = SO3(dR).inverse() * (R1*SO3::exp(dd_phi)).inverse() * R2;
    //数值方法计算这一步的导数
    Mat3d result = (eR_dd_phi - eR)/dd_phi;
    //对比result和 Mat3d -invJr * (R2.inverse() * R1).matrix()

    //自己修改的部分end
    
    //旋转部分残差
    const Vec3d er = eR.log();
    const Mat3d invJr = SO3::jr_inv(eR);

_jacobianOplus[0].block<3, 3>(0, 0) = -invJr * (R2.inverse() * R1).matrix();

//自己写的数值求导方法
//根据4.41a的定义

