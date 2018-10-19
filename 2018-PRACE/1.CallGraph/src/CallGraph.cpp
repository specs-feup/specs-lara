namespace test {

	class B {
	};

	class A {
		public:
		B foo();
	};

}

void callFoo(test::B);

int main() {
	test::A a;
	
	callFoo(a.foo());
}