import geometry.Circle;
import geometry.Rectangle;
import screen.Screen;

public class ApplicationObserver {
	public static void main(String[] args) {
		Screen screen = new Screen();
		
		Rectangle rectangle = new Rectangle(10, 10 ,20 ,20);
		Circle circle = new Circle(15, 15 , 10);
		
		screen.addShape(rectangle);
		screen.addShape(circle);
		
		rectangle.setX(rectangle.getX() + 10);
		
		screen.removeShape(circle);
	}
}