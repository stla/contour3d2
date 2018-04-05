module ExeOrthocircle
  where
import           Lib
import           Text.Show.Pretty

f :: Double -> Double -> Double -> Double
f x y z =
    ((x2 + y2 - 1)**2 + z2) * ((y2 + z2 - 1)**2 + x2) *
        ((z2 + x2 - 1)**2 + y2) - a**2*(1 + b*(x2 + y2 + z2))
    where
    a = 0.075
    b = 3
    x2 = x*x
    y2 = y*y
    z2 = z*z

triangles :: [[Double]]
triangles = computeContour3d (fun2array 300 (-1.6,1.6) f) Nothing 0 Nothing Nothing Nothing

main :: IO ()
main = writeFile "zorthocircle300_show.txt" (show triangles)
