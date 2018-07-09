package screen;

import java.util.LinkedList;
import java.util.List;

import geometry.Shape;

public class Screen {
	List<Shape> shapes = new LinkedList<>();
	
	public void addShape(Shape shape) {
		shapes.add(shape);
		draw();
	}

	public  void removeShape(Shape shape) {
		shapes.remove(shape);
		draw();
	}
	
	public void draw() {
		// Method that would render the shapes to the screen
		System.out.println("Drawing " + shapes.size() + " shapes");
	}
}
