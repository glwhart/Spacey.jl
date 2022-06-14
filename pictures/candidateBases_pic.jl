#using Luxor
begin
Drawing(500, 500, "my-drawing.svg")
origin()
setcolor("red")
circle(Point(0, 0), 100, :fill)
finish()
end

begin
    Drawing(200, 200, "my-drawing.png")
background(0.9, 0.2, 0.2, .5)
setcolor("blue")
setline(10)
setdash("dash")
circle(Point(0, 0), 100, :stroke)
finish()

end