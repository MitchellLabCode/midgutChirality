function myKeyPressCallback(~, event)
    global restartLoopFlag;
    global finishedFlag;
    if strcmp(event.Key, 'r') % Check if 'r' key was pressed
        restartLoopFlag = true;
        disp('R key pressed. Setting restart flag.');
    end
    if strcmp(event.Key, 'f') % Check if 'f' key was pressed
        finishedFlag = true;
        disp('f key pressed. Next frame.');
    end
end