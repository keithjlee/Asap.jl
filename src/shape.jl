```
SHAPE FUNCTIONS gives displacement orthogonal to local x axis based on end moments and shear forces
u1 = shear displacement at beginning node
u2 = moment displacement at beginning node
u3 = shear displacement at end node
u4 = moment displacement at end node

These displacements should be w/r/t local coordinate system
```

# 
function N1(x, element::Element)
    return 1 - 3 * (x / element.length)^2 + 2 * (x / element.length)^3
end

function N2(x, element::Element)
    return x * (1 - x / element.length)^2
end

function N3(x, element::Element)
    return 3 * (x / element.length)^2 - 2 * (x / element.length)^3
end

function N4(x, element::Element)
    return x^2 / element.length * (x / element.length - 1)
end