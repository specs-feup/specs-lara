import geometry.Circle;
import geometry.Rectangle;

public class ApplicationLogging {
	public static void main(String[] args) {
		Rectangle rectangle = new Rectangle(10, 10 ,20 ,20);
		Circle circle = new Circle(15, 15 , 10);
		
		rectangle.setX(rectangle.getX() + 10);
		circle.setRadius(20);
	}
}