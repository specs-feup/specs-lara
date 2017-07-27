function dynamicCallGraph
	global call_graph
	
	call_graph
	
	fprintf('digraph CallGraph {\n');
		if (numel(call_graph) >= 1 && call_graph(1) ~= 0) 
			fprintf('\tfoobar -> foo [label=\"%d\"];\n', call_graph(1));
		end
		if (numel(call_graph) >= 2 && call_graph(2) ~= 0)
			fprintf('\tfoobar -> other [label=\"%d\"];\n', call_graph(2));
		end
		if (numel(call_graph) >= 3 && call_graph(3) ~= 0)
			fprintf('\tfoo -> bar [label=\"%d\"];\n', call_graph(3));
		end

	fprintf('}\n');

	call_graph = [];
end