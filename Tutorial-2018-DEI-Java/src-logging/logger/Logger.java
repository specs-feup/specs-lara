package logger;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Logger {
	private static String INFO = "I";
	private static String WARNING = "W";
	private static String DEBUG = "D";
	private static String ERROR = "E";
	
	private static void log(String level, String message) {
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		
		System.out.println("| " + level + " | " + dateFormat.format(date) + " | " + message + " |");
	}
	
	public static void i(String message) {
		log(INFO, message);
	}

	public static void w(String message) {
		log(WARNING, message);
	}

	public static void d(String message) {
		log(DEBUG, message);
	}

	public static void e(String message) {
		log(ERROR, message);
	}
}
