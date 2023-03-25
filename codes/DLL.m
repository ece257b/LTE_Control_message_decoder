classdef DLL < handle
    % DLL an element of a doubly-linked list.
    
    properties
        val  % data stored in this element
    end
    
    properties (Access = private)
        next % reference to next DLL object
        prev % reference to previous DLL object
    end
    
    methods
        function obj = DLL(val)
        obj.val = val;
        obj.next = obj;
        obj.prev = obj;
        end
        
        function obj = pop(obj)
        obj.prev.next = obj.next;
        obj.next.prev = obj.prev;
        obj.next = obj;
        obj.prev = obj;
        end
        
        function obj = append(obj, other)
        other.next = obj.next;
        other.prev = obj;
        obj.next.prev = other;
        obj.next = other;
        end
        
        function other = get_next(obj)
        other = obj.next;
        end
        
        function other = get_prev(obj)
        other = obj.prev;
        end
    end
end