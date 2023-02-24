function dcm_from_mrp(p)
    p1,p2,p3 = p
    den = (p1^2 + p2^2 + p3^2 + 1)^2
    a = (4*p1^2 + 4*p2^2 + 4*p3^2 - 4)
    [
    (-((8*p2^2+8*p3^2)/den-1)*den)   (8*p1*p2 + p3*a)     (8*p1*p3 - p2*a);
    (8*p1*p2 - p3*a) (-((8*p1^2 + 8*p3^2)/den - 1)*den)   (8*p2*p3 + p1*a);
    (8*p1*p3 + p2*a)  (8*p2*p3 - p1*a)  (-((8*p1^2 + 8*p2^2)/den - 1)*den)
    ]/den
end
function skew(ω::Vector{T}) where {T}
    return [0 -ω[3] ω[2];
            ω[3] 0 -ω[1];
            -ω[2] ω[1] 0]
end
function dynamics(model::NamedTuple,x,u)
    # quadrotor dynamics with an MRP for attitude
    r = x[1:3]
    v = x[4:6]
    p = x[7:9]
    ω = x[10:12]

    Q = dcm_from_mrp(p)

    mass=model.mass
    J = model.J
    gravity= model.gravity
    L= model.L
    kf=model.kf
    km=model.km

    w1 = u[1]
    w2 = u[2]
    w3 = u[3]
    w4 = u[4]

    F1 = max(0,kf*w1)
    F2 = max(0,kf*w2)
    F3 = max(0,kf*w3)
    F4 = max(0,kf*w4)
    F = [0., 0., F1+F2+F3+F4] #total rotor force in body frame

    M1 = km*w1
    M2 = km*w2
    M3 = km*w3
    M4 = km*w4
    τ = [L*(F2-F4), L*(F3-F1), (M1-M2+M3-M4)] #total rotor torque in body frame

    f = mass*gravity + Q*F # forces in world frame

    [
        v
        f/mass
        ((1+norm(p)^2)/4) *(   I + 2*(skew(p)^2 + skew(p))/(1+norm(p)^2)   )*ω
        J\(τ - cross(ω,J*ω))
    ]
end
function rk4(model,ode,x,u,dt)
    k1 = dt*ode(model,x, u)
    k2 = dt*ode(model,x + k1/2, u)
    k3 = dt*ode(model,x + k2/2, u)
    k4 = dt*ode(model,x + k3, u)
    x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end
function vis_traj!(vis, name, X; R = 0.1, color = mc.RGBA(1.0, 0.0, 0.0, 1.0))
    for i = 1:(length(X)-1)
        a = X[i][1:3]
        b = X[i+1][1:3]
        cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
        mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color))
    end
    for i = 1:length(X)
        a = X[i][1:3]
        sph = mc.HyperSphere(mc.Point(a...), R)
        mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color))
    end
end

function animate_quadrotor(Xsim, Xref, dt)
    vis = mc.Visualizer()
    robot_obj = mc.MeshFileGeometry(joinpath(@__DIR__,"quadrotor.obj"))
    mc.setobject!(vis[:vic], robot_obj)

    vis_traj!(vis, :traj, Xref[1:85]; R = 0.01, color = mc.RGBA(1.0, 0.0, 0.0, 1.0))
    target = mc.HyperSphere(mc.Point(0,0,0.0),0.1)
    mc.setobject!(vis[:target], target, mc.MeshPhongMaterial(color = mc.RGBA(0.0,1.0,0.0,0.4)))


    anim = mc.Animation(floor(Int,1/dt))
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            r = Xsim[k][1:3]
            p = Xsim[k][7:9]
            mc.settransform!(vis[:vic], mc.compose(mc.Translation(r),mc.LinearMap(1.5*(dcm_from_mrp(p)))))
            mc.settransform!(vis[:target], mc.Translation(Xref[k][1:3]))
        end
    end
    mc.setanimation!(vis, anim)

    return (mc.render(vis))
end