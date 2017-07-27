function foobar()
    global call_graph
    
	foo();
	if(numel(call_graph) < 1); call_graph(1) = 1; else; call_graph(1) = call_graph(1) + 1; end
	other();
	if(numel(call_graph) < 2); call_graph(2) = 1; else; call_graph(2) = call_graph(2) + 1; end
end

function a = foo()
	global call_graph
    
	a = 0;
    
    for i = 1:5
        a = a + bar();
		if(numel(call_graph) < 3); call_graph(3) = 1; else; call_graph(3) = call_graph(3) + 1; end
    end
    
end


function result = bar() 
	global call_graph
    
	result = 1.0;
end


