package person;

import java.util.LinkedList;
import java.util.List;

public class Person {
    /***
     * A list of persons in the system.
     */
    private static List<Person> people = new LinkedList<>();

    /***
     * Add a new person to the system.
     */
    public static void addPerson(Person person) {
        people.add(person);
    }

    /**
     * Get all the people in the system.
     * 
     * @return All the people in the system.
     */
    public static List<Person> getPeople() {
        return people;
    }

    /**
     * The name of this person.
     */
    private String name;

    /**
     * The age of this person.
     */
    private int age;

    /**
     * A constructor taking the name and age of a person.
     */
    public Person(String name, int age) {
        this.name = name;
        this.age = age;
    }

    /**
     * Returns the name of this person.
     * 
     * @return This person's name.
     */
    public String getName() {
        return name;
    }

    /**
     * Changes the name of this person.
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Returns the age of this person.
     * 
     * @return This person's age.
     */
    public int getAge() {
        return age;
    }

    /**
     * Changes the age of this person.
     */
    public void setAge(int age) {
        this.age = age;
    }
}