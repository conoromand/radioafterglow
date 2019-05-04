references：
Granot, Piran, and Sari 99 (GPS99)
Granot and Piran 12 (GP12)

大まかな計算の流れ：
1. GP12のモデルを使ってjetのdynamicsを計算 --> Gamma(t), R(t), theta(t)
2. 1.からshocked regionの電子加速とsync.emissivityを計算 --> P'_nu(nu',r,t) (注: 各tでtheta(t)の外側ではP'_nu(nu',r,t) = 0)
3. GPS99のEq. (2)に従って、given T_obsに対する"egg shape"を求める。
4. 3.で求めたegg shapeに従ってGPS99のEq. (4)の積分を実行。
5. 4.を繰り返してlight curve or spectrumを計算。
