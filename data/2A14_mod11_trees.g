Tree_3 := rec(
    name := "[3]",
    center := [
        A14SpinPos(t, [11,3], '+'),
        A14SpinPos(t, [11,3], '-')
    ],
    branches := [
        [
            A14SpinPos(t, [10,3,1], 0),
            A14SpinPos(t, [9,3,2], 0),
            A14SpinPos(t, [7,4,3], 0),
            A14SpinPos(t, [6,5,3], 0)
        ],
        [
            A14SpinPos(t, [14], 0)
        ]
    ],
    exceptionalMultiplicity := 2
);

Tree_21 := rec(
    name := "[2,1]",
    center := [
        A14SpinPos(t, [11,2,1], 0)
    ],
    branches := [
        [
            A14SpinPos(t, [12,2], '+'),
            A14SpinPos(t, [13,1], '+')
        ],
        [
            A14SpinPos(t, [8,3,2,1], '+'),
            A14SpinPos(t, [7,4,2,1], '+'),
            A14SpinPos(t, [6,5,2,1], '+')
        ],
        [
            A14SpinPos(t, [12,2], '-'),
            A14SpinPos(t, [13,1], '-')
        ],
        [
            A14SpinPos(t, [8,3,2,1], '-'),
            A14SpinPos(t, [7,4,2,1], '-'),
            A14SpinPos(t, [6,5,2,1], '-')
        ]
    ],
    exceptionalMultiplicity := 1
);

trees := [ Tree_3, Tree_21 ];