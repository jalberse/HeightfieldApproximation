# Heightfield Approximation

![gif of a build](media/interact.gif)

An interactive visualization of the heightfield approximation algorithm for realtime fluid surface simulation over arbitrary domains, built on top of the [PixelGameEngine](https://github.com/OneLoneCoder/olcPixelGameEngine/wiki/olc::PixelGameEngine).

Heightfield approximation is a computationally cheap and relatively simple method which gives surprisingly realistic looking results. In exchange, the approximation is unable to represent more complex phenomena such as breaking waves, splashing, and droplets. This makes the algorithm a good fit for real-time simulation of large bodies of fluid such as lakes or oceans. 

Currently, this visualization represents the heightfield in 2D, however it is possible to use the data to generate a polygonal mesh to view the fluid surface in 3D. 

| Option | Key |
|---|---|
| Pause/Resume |  Space |
| Perturb fluid surface | Left Click  |
| Add terrain cell | Right Click  |
| Remove terrain tell | Shift + Right Click |
| Clear terrain | C |
| Decrease/Increase dampening term | N/M |
| Reset fluid to new perlin noise | R |
| Reset fluid and terrain to new perlin noise | Shift + R |
| Decrease/Increase perlin octave | O/P |
| Decrease/Increase perlin scaling Bias | K/L | 
| Decrease/Increase perlin octave (terrain) | Shift + O/P |
| Decrease/Increase perlin scaling bias (Terrain) | Shift + K/L |