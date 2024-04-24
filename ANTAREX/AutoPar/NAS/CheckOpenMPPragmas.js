laraImport("lara.Io");
laraImport("lara.Check");
laraImport("weaver.Query");

function CheckOpenMPPragmas(expectedCodeFile, generateOutputs = false) {
  println("Checking OpenMP pragmas...");
  //	select pragma{"omp"}.target end
  //	apply
  for (const $pragma of Query.search("pragma", "omp")) {
    const $target = $pragma.target;

    // Comment OpenMP pragmas that are inside a loop that already has an OpenMP pragma
    if ($target.instanceOf("loop")) {
      if (!hasForAncestorWithOmp($target)) {
        continue;
      }

      $pragma.insertReplace("// " + $pragma.code);
    }
  }
  //	end

  //select file end
  //apply
  let currentCode = "<not initialized>";
  for (const $file of Query.search("file")) {
    currentCode = $file.code;
    break;
  }
  //	end

  // Generate outputs and return
  if (generateOutputs) {
    println("Generating outputs and returning");
    Io.writeFile(expectedCodeFile, currentCode);
    return;
  }

  if (expectedCodeFile === undefined) {
    println("No expected code file, returning");
    return;
  }

  var expectedCode = Io.readFile(expectedCodeFile);
  if (expectedCode === null) {
    expectedCode = "Could not find file '" + expectedCodeFile + "'";
  }

  Check.strings(currentCode, expectedCode);
  println("Code is as expected");
}

function hasForAncestorWithOmp($target) {
  // Find ancestor that is a loop
  $loopParent = $target.ancestor("loop");

  // No loop parent found, return
  if ($loopParent === undefined) {
    return false;
  }

  // Check if parent has an OpenMP pragma
  for (var $parentPragma of $loopParent.pragmas) {
    // Found OpenMP pragma, return
    if ($parentPragma.name === "omp") {
      return true;
    }
  }

  // No OpenMP pragma found, check parent for
  return hasForAncestorWithOmp($loopParent);
}

const expectedCodeFile = laraArgs["expectedCodeFile"];
const generateOutputs = laraArgs["generateOutputs"];
CheckOpenMPPragmas(expectedCodeFile, generateOutputs);
