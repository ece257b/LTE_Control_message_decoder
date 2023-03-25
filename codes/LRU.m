classdef LRU < handle
    %LRU a least-recently-used cache.
    
    properties (Access = private)
        sentinel  % sentinel of doubly-linked-list containing cached data. sentinel.next is most-
        % recently used, sentinel.prev is least-recently used
        map       % map of key --> dll so that dll.val is the value associated with the given key
        itemlimit % maximum size of cache in number of items
        memlimit  % maximum size of cache in bytes (default inf)
        mem       % current memory usage in all cached objects (bytes)
    end
    
    methods
        function obj = LRU(itemlimit, memlimit)
        if nargin < 2, memlimit = inf; end
        
        % Initialize DLL
        obj.sentinel = DLL([]);
        
        % Initialize Map
        obj.map = containers.Map();
        
        % Set limits
        obj.itemlimit = itemlimit;
        obj.memlimit = memlimit;
        
        % Initial memory usage is 0
        obj.mem = 0;
        end
        
        function val = get(obj, key)
        if obj.iscached(key)
            dll = obj.map(key);
            val = dll.val{2};
            
            % Re-insert dll at the front of the list since it is now most recently used
            dll.pop();
            obj.sentinel.append(dll);
        else
            val = [];
        end
        end
        
        function val = put(obj, key, val)
        % If key is already in map, remove existing object (need to update memory count, etc)
        if obj.iscached(key)
            obj.remove(key);
        end
        
        % Wrap value in a DLL node and insert into the 'most recently used' position.
        dll = DLL({key, val});
        obj.sentinel.append(dll);
        obj.map(key) = dll;
        
        % Increment total memory usage count
        obj.mem = obj.mem + bytes(val);
        
        % Ensure that memory is kept under limit by dropping old items until we are back below the
        % limit (memory or total items limits)
        while obj.mem > obj.memlimit || length(obj.map) > obj.itemlimit
            obj.drop_oldest();
        end
        end
        
        function val = remove(obj, key)
        if obj.iscached(key)
            % Remove object from DLL
            dll = obj.map(key);
            dll.pop();
            
            % Decrement memory count
            obj.mem = obj.mem - bytes(dll.val{2});
            
            % Remove key from map
            obj.map.remove(key);
        else
            % Nothing to do
            val = [];
        end
        end
        
        function tf = iscached(obj, key)
            tf = obj.map.isKey(key);
        end
        
        function s = size(obj)
        s = length(obj.map);
        end
        
        function m = memusage(obj)
        m = obj.mem;
        end
    end
    
    methods (Access=private)
        function drop_oldest(obj)
        oldest = obj.sentinel.get_prev();
        oldest.pop();
        obj.mem = obj.mem - bytes(oldest.val{2});
        obj.map.remove(oldest.val{1});
        end
    end
    
end

function b = bytes(var)
info = whos('var');
b = info.bytes;
end