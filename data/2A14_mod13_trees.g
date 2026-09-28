Tree_1 := rec(
    name := "[1]",
    center := [
        A14SpinPos(t, [13,1], '+'),
        A14SpinPos(t, [13,1], '-')
    ],
    branches := [
        [
            A14SpinPos(t, [7,6,1], 0),
            A14SpinPos(t, [8,5,1], 0),
            A14SpinPos(t, [9,4,1], 0),
            A14SpinPos(t, [10,3,1], 0),
            A14SpinPos(t, [11,2,1], 0)
        ],
        [
            A14SpinPos(t, [14], 0)
        ]
    ],
    exceptionalMultiplicity := 2
);

trees := [ Tree_1 ];