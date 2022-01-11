# RB_rdyn_4unity

Implementation of exact rotation dynamics for Unity

## Unity Physix is wrong

When it comes to simulating rotations there are some unintuitiv behaviors in the real world. 
[Youtube Video from Veritasium](https://youtu.be/1VPfZ_XzisU)

There are 2 things to consider for an object floating undisturbed in free space: 
* Conservation of angular momentum and 
* Conservation of kinetic energy.

But these do not state that the "internal" orientation of the body must stay the same while it rotates in the momentum axis. 
This goes further that the accelerations of the bodys mass distributions might even demand a change of the "internal" rotation.

This effect highly depends on the mass distribution - more commonly known its **inertia**.
  
![before flip](/doc/flip1.png)
![after flip](/doc/flip2.png)

A standard Unity simulation with Physix will progress the rotation of an object just as a continuos 
rotation around the fixed axis. This is correct if the principal inertias are all equal. As soon as 
the inertias are non-equal this is not correct any more.  

This implementation demonstrates the non-intuitive rotations of arbitary inertia distributions.
This algorithm provides an **analytic solution**. It does not need to advance regularly - 
you can jump years ahead and get the correct rotation at this exact moment. Arbitary time movements can be done.

The algorithm was implemented with reference to the paper [1][1].
The underlying code uses the jacqobi elliptic functions to look up the rotation axis for a specific time t.

    ellipk.ellpj(u, m, ref sn, ref cn, ref dn, ref ph);

    omegas.x = wm_s.x * cn;
    omegas.y = wm_s.y * sn;
    omegas.z = wm_s.z * dn;

Elliptic functions were translated from the cephes library [2][2] used also as the scipy reference.


## How to use

On a gameobject attach the script `source/FreeRotor.cs`

Set the principal inertia components to something > 0. 
Set the angular velocity components to something > 0.



Maybe try to avoid unrealistic inertia settings. 
For a rod-like shape you could use (1,10,10.1).
For a quad something like (5,7,7.5), or (5,5.5,9).


### Limitations

Currently only the setting of principal inertia is allowed. This means your object in body space should be symmetric around each individual axis.
I think there is some bug with specific ang-velocity settings, it seems to accelerate and brake the momentum.
Kin-Energy conservation results not always static - maybe bug somewhere?

## Todos

What I want to add:
* approximate inertia matrix from mesh
* initial rotation to form pricipal inertias
* add cases for symmetric and spherical tops to skip complicated movements
* default inertias if set to 0
* testing for very high / very low numbers
* check momentum conservation

Ideas:
* some others


## References

[1]: https://doi.org/10.1016/j.jcp.2006.11.019 (Numerical implementation of the exact dynamics of free rigid bodies)  
[2]: http://www.netlib.org/cephes/ (Cephes function library)  

## License

**[MIT](LICENSE)**