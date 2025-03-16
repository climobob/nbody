# nbody
Simple symplectic integrator for the N body problem

solar_system.py includes part of the solar system and no moons -- Sun, Venus, Earth, Jupiter, Saturn
run by python3 solar_system.py > output
and get a text file with the locations (x,y) of each body and its change in distance from initial distance to the system center of mass

earth_system.py is for an earth with a number of lunar mass bodies orbiting it. One is at 1 lunar distance, at earth-moon Keplerian speed. The rest are uniformly random distances in range [0.3,3] with velocities uniformly random between [0.2, sqrt(2) ) times the Keplerian speed. The sun is ignored in this calculation.
run by python3 earth_system.py > output
same output format as solar_system.py.  The number of moons is not constant. They are allowed to collide and merge.





